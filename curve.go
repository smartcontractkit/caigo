package caigo

/*
	Although the library adheres to the 'elliptic/curve' interface.
	All testing has been done against library function explicity.
	It is recommended to use in the same way(i.e. `curve.Sign` and not `ecdsa.Sign`).
*/
import (
	"crypto/elliptic"
	_ "embed"
	"encoding/json"
	"fmt"
	"log"
	"math/big"
)

var Curve StarkCurve

/*
Returned stark curve includes several values above and beyond
what the 'elliptic' interface calls for to facilitate common starkware functions
*/
type StarkCurve struct {
	*elliptic.CurveParams
	EcGenX           *big.Int
	EcGenY           *big.Int
	MinusShiftPointX *big.Int
	MinusShiftPointY *big.Int
	Max              *big.Int
	Alpha            *big.Int
	ConstantPoints   [][]*big.Int
}

//go:embed pedersen_params.json
var PedersenParamsRaw []byte
var PedersenParams StarkCurvePayload

// struct definition for parsing 'pedersen_params.json'
type StarkCurvePayload struct {
	License        []string     `json:"_license"`
	Comment        string       `json:"_comment"`
	FieldPrime     *big.Int     `json:"FIELD_PRIME"`
	FieldGen       int          `json:"FIELD_GEN"`
	EcOrder        *big.Int     `json:"EC_ORDER"`
	Alpha          int64        `json:"ALPHA"`
	Beta           *big.Int     `json:"BETA"`
	ConstantPoints [][]*big.Int `json:"CONSTANT_POINTS"`
}

func init() {
	if err := json.Unmarshal(PedersenParamsRaw, &PedersenParams); err != nil {
		log.Fatalf("unmarshalling pedersen params: %v", err)
	}

	if len(PedersenParams.ConstantPoints) == 0 {
		panic("decoding pedersen params json")
	}

	Curve.CurveParams = &elliptic.CurveParams{Name: "stark-curve-with-constants"}

	Curve.P = PedersenParams.FieldPrime
	Curve.N = PedersenParams.EcOrder
	Curve.B = PedersenParams.Beta
	Curve.Gx = PedersenParams.ConstantPoints[0][0]
	Curve.Gy = PedersenParams.ConstantPoints[0][1]
	Curve.EcGenX = PedersenParams.ConstantPoints[1][0]
	Curve.EcGenY = PedersenParams.ConstantPoints[1][1]
	Curve.MinusShiftPointX, _ = new(big.Int).SetString("2089986280348253421170679821480865132823066470938446095505822317253594081284", 10) // MINUS_SHIFT_POINT = (SHIFT_POINT[0], FIELD_PRIME - SHIFT_POINT[1])
	Curve.MinusShiftPointY, _ = new(big.Int).SetString("1904571459125470836673916673895659690812401348070794621786009710606664325495", 10)
	Curve.Max, _ = new(big.Int).SetString("3618502788666131106986593281521497120414687020801267626233049500247285301248", 10) // 2 ** 251
	Curve.Alpha = big.NewInt(PedersenParams.Alpha)
	Curve.BitSize = 252
	Curve.ConstantPoints = PedersenParams.ConstantPoints

	/*
		Not all operations require a stark curve initialization
		including the provided constant points. Here you can
		initialize the curve without the constant points
	*/
	Curve.CurveParams = &elliptic.CurveParams{Name: "stark-curve"}
	Curve.P, _ = new(big.Int).SetString("3618502788666131213697322783095070105623107215331596699973092056135872020481", 10)  // Field Prime ./pedersen_json
	Curve.N, _ = new(big.Int).SetString("3618502788666131213697322783095070105526743751716087489154079457884512865583", 10)  // Order of base point ./pedersen_json
	Curve.B, _ = new(big.Int).SetString("3141592653589793238462643383279502884197169399375105820974944592307816406665", 10)  // Constant of curve equation ./pedersen_json
	Curve.Gx, _ = new(big.Int).SetString("2089986280348253421170679821480865132823066470938446095505822317253594081284", 10) // (x, _) of basepoint ./pedersen_json
	Curve.Gy, _ = new(big.Int).SetString("1713931329540660377023406109199410414810705867260802078187082345529207694986", 10) // (_, y) of basepoint ./pedersen_json
	Curve.EcGenX, _ = new(big.Int).SetString("874739451078007766457464989774322083649278607533249481151382481072868806602", 10)
	Curve.EcGenY, _ = new(big.Int).SetString("152666792071518830868575557812948353041420400780739481342941381225525861407", 10)
	Curve.MinusShiftPointX, _ = new(big.Int).SetString("2089986280348253421170679821480865132823066470938446095505822317253594081284", 10) // MINUS_SHIFT_POINT = (SHIFT_POINT[0], FIELD_PRIME - SHIFT_POINT[1])
	Curve.MinusShiftPointY, _ = new(big.Int).SetString("1904571459125470836673916673895659690812401348070794621786009710606664325495", 10) // MINUS_SHIFT_POINT = (SHIFT_POINT[0], FIELD_PRIME - SHIFT_POINT[1])
	Curve.Max, _ = new(big.Int).SetString("3618502788666131106986593281521497120414687020801267626233049500247285301248", 10)              // 2 ** 251
	Curve.Alpha = big.NewInt(1)
	Curve.BitSize = 252
}

// Gets two points on an elliptic curve mod p and returns their sum.
// Assumes affine form (x, y) is spread (x1 *big.Int, y1 *big.Int)
//
// (ref: https://github.com/starkware-libs/cairo-lang/blob/master/src/starkware/crypto/starkware/crypto/signature/math_utils.py)
func (sc StarkCurve) Add(x1, y1, x2, y2 *big.Int) (x, y *big.Int) {
	// As elliptic curves form a group, there is an additive identity that is the equivalent of 0
	// If 𝑃=0 or 𝑄=0, then 𝑃+𝑄=𝑄 or 𝑃+𝑄=𝑃, respectively
	// NOTICE: the EC multiplication algorithm is using using `StarkCurve.rewriteScalar` trick
	//   to avoid this condition and provide constant-time execution.
	if len(x1.Bits()) == 0 && len(y1.Bits()) == 0 {
		return x2, y2
	}
	if len(x2.Bits()) == 0 && len(y2.Bits()) == 0 {
		return x1, y1
	}

	yDelta := new(big.Int).Sub(y1, y2)
	xDelta := new(big.Int).Sub(x1, x2)

	m := DivMod(yDelta, xDelta, sc.P)

	xm := new(big.Int).Mul(m, m)

	x = new(big.Int).Sub(xm, x1)
	x = x.Sub(x, x2)
	x = x.Mod(x, sc.P)

	y = new(big.Int).Sub(x1, x)
	y = y.Mul(m, y)
	y = y.Sub(y, y1)
	y = y.Mod(y, sc.P)

	return x, y
}

// Doubles a point on an elliptic curve with the equation y^2 = x^3 + alpha*x + beta mod p.
// Assumes affine form (x, y) is spread (x1 *big.Int, y1 *big.Int)
//
// (ref: https://github.com/starkware-libs/cairo-lang/blob/master/src/starkware/crypto/starkware/crypto/signature/math_utils.py)
func (sc StarkCurve) Double(x1, y1 *big.Int) (x, y *big.Int) {
	xin := new(big.Int).Mul(big.NewInt(3), x1)
	xin = xin.Mul(xin, x1)
	xin = xin.Add(xin, sc.Alpha)

	yin := new(big.Int).Mul(y1, big.NewInt(2))

	m := DivMod(xin, yin, sc.P)

	xout := new(big.Int).Mul(m, m)
	xmed := new(big.Int).Mul(big.NewInt(2), x1)
	xout = xout.Sub(xout, xmed)
	xout = xout.Mod(xout, sc.P)

	yout := new(big.Int).Sub(x1, xout)
	yout = yout.Mul(m, yout)
	yout = yout.Sub(yout, y1)
	yout = yout.Mod(yout, sc.P)

	return xout, yout
}

func (sc StarkCurve) ScalarMult(x1, y1 *big.Int, k []byte) (x, y *big.Int) {
	m := new(big.Int).SetBytes(k)
	x, y = sc.EcMult(m, x1, y1)
	return x, y
}

func (sc StarkCurve) ScalarBaseMult(k []byte) (x, y *big.Int) {
	return sc.ScalarMult(sc.Gx, sc.Gy, k)
}

func (sc StarkCurve) IsOnCurve(x, y *big.Int) bool {
	left := new(big.Int).Mul(y, y)
	left = left.Mod(left, sc.P)

	right := new(big.Int).Mul(x, x)
	right = right.Mul(right, x)
	right = right.Mod(right, sc.P)

	ri := new(big.Int).Mul(big.NewInt(1), x)

	right = right.Add(right, ri)
	right = right.Add(right, sc.B)
	right = right.Mod(right, sc.P)

	return left.Cmp(right) == 0
}

// (ref: https://github.com/starkware-libs/cairo-lang/blob/master/src/starkware/crypto/starkware/crypto/signature/math_utils.py)
func (sc StarkCurve) InvModCurveSize(x *big.Int) *big.Int {
	return DivMod(big.NewInt(1), x, sc.N)
}

// Given the x coordinate of a stark_key, returns a possible y coordinate such that together the
// point (x,y) is on the curve.
// Note: the real y coordinate is either y or -y.
//
// (ref: https://github.com/starkware-libs/cairo-lang/blob/master/src/starkware/crypto/starkware/crypto/signature/signature.py)
func (sc StarkCurve) GetYCoordinate(starkX *big.Int) *big.Int {
	y := new(big.Int).Mul(starkX, starkX)
	y = y.Mul(y, starkX)
	yin := new(big.Int).Mul(sc.Alpha, starkX)

	y = y.Add(y, yin)
	y = y.Add(y, sc.B)
	y = y.Mod(y, sc.P)

	y = y.ModSqrt(y, sc.P)
	return y
}

// Computes m * point + shift_point using the same steps like the AIR and throws an exception if
// and only if the AIR errors.
//
// (ref: https://github.com/starkware-libs/cairo-lang/blob/master/src/starkware/crypto/starkware/crypto/signature/signature.py)
func (sc StarkCurve) MimicEcMultAir(mout, x1, y1, x2, y2 *big.Int) (x *big.Int, y *big.Int, err error) {
	m := new(big.Int).Set(mout)
	if m.Cmp(big.NewInt(0)) != 1 || m.Cmp(sc.Max) != -1 {
		return x, y, fmt.Errorf("too many bits %v", m.BitLen())
	}

	psx := x2
	psy := y2
	for i := 0; i < 251; i++ {
		if psx == x1 {
			return x, y, fmt.Errorf("xs are the same")
		}
		if m.Bit(0) == 1 {
			psx, psy = sc.Add(psx, psy, x1, y1)
		}
		x1, y1 = sc.Double(x1, y1)
		m = m.Rsh(m, 1)
	}
	if m.Cmp(big.NewInt(0)) != 0 {
		return psx, psy, fmt.Errorf("m doesn't equal zero")
	}
	return psx, psy, nil
}

// Multiplies by m a point on the elliptic curve with equation y^2 = x^3 + alpha*x + beta mod p.
// Assumes affine form (x, y) is spread (x1 *big.Int, y1 *big.Int) and that 0 < m < order(point).
//
// (ref: https://www.semanticscholar.org/paper/Elliptic-Curves-and-Side-Channel-Analysis-Joye/7fc91d3684f1ab63b97d125161daf57af60f2ad9/figure/1)
// (ref: https://cosade.telecom-paristech.fr/presentations/s2_p2.pdf)
func (sc StarkCurve) ecMult_DoubleAndAlwaysAdd(m, x1, y1 *big.Int) (x, y *big.Int) {
	var _ecMult = func(m, x1, y1 *big.Int) (x, y *big.Int) {
		// Two-index table initialization, Q[0], Q[1] <- P, 2*P
		dx1, dy1 := sc.Double(x1, y1)
		q := [2]struct {
			x *big.Int
			y *big.Int
		}{
			{
				x: x1,
				y: y1,
			},
			{
				x: dx1,
				y: dy1,
			},
		}

		// Run the algorithm, expects the most-significant bit is 1 (skip MSB)
		for i := sc.BitSize - 1; i >= 0; i-- {
			b := m.Bit(i)                                               // b      <- bit at position i
			q[1-b].x, q[1-b].y = sc.Add(q[0].x, q[0].y, q[1].x, q[1].y) // Q[1-b] <- Q[0] + Q[1]
			q[b].x, q[b].y = sc.Double(q[b].x, q[b].y)                  // Q[b]   <- 2Q[b]
		}

		return q[0].x, q[0].y
	}

	return _ecMult(sc.rewriteScalar(m), x1, y1)
}

// Rewrites k into an equivalent scalar, such that the first bit (the most-significant
// bit for the Double-And-Always-Add or Montgomery algo) is 1.
//
// The k scalar rewriting obtains an equivalent scalar K = 2^n + (k - 2^n mod q),
// such that k·G == K·G and K has the n-th bit set to 1. The scalars are equal modulo
// the group order, k mod q == K mod q.
//
// Notice: The EC multiplication algorithms are typically presented as starting with the state (O, P0),
//   where O is the identity element (or neutral point) of the curve. However, the neutral point is at infinity,
//   which causes problems for some formulas (non constant-time execution for the naive implementation).
//   The ladder then starts after the first step, when the state no longer contains the neutral point.
// (ref: https://www.shiftleft.org/papers/ladder/ladder-tches.pdf)
func (sc StarkCurve) rewriteScalar(k *big.Int) *big.Int {
	size := new(big.Int).Lsh(big.NewInt(1), uint(sc.BitSize)) // 2ˆn
	diff := new(big.Int).Sub(k, size)                         // (k - 2ˆn)
	diffMod := new(big.Int).Mod(diff, sc.N)                   // (k - 2ˆn) mod q
	return new(big.Int).Add(size, diffMod)                    // 2ˆn + (k - 2ˆn mod q)
}

// Multiplies by m a point on the elliptic curve with equation y^2 = x^3 + alpha*x + beta mod p.
// Assumes affine form (x, y) is spread (x1 *big.Int, y1 *big.Int) and that 0 < m < order(point).
func (sc StarkCurve) EcMult(m, x1, y1 *big.Int) (x, y *big.Int) {
	return sc.ecMult_DoubleAndAlwaysAdd(m, x1, y1)
}

// Finds a nonnegative integer 0 <= x < p such that (m * x) % p == n
//
// (ref: https://github.com/starkware-libs/cairo-lang/blob/master/src/starkware/crypto/starkware/crypto/signature/math_utils.py)
func DivMod(n, m, p *big.Int) *big.Int {
	q := new(big.Int)
	gx := new(big.Int)
	gy := new(big.Int)
	q.GCD(gx, gy, m, p)

	r := new(big.Int).Mul(n, gx)
	r = r.Mod(r, p)
	return r
}
