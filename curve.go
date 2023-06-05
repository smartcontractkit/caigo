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
	B3               *big.Int
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

	Curve.B3 = new(big.Int).Mul(Curve.B, big.NewInt(3))
	Curve.B3.Mod(Curve.B3, Curve.P)
}

// Gets two points on an elliptic curve mod p and returns their sum.
// Assumes affine form (x, y) is spread (x1 *big.Int, y1 *big.Int)
//
// This implements Algorithm 2 of 2015 Renes–Costello–Batina "Complete addition formulas for prime order elliptic curves"
// (ref: https://eprint.iacr.org/2015/1060.pdf)
func (sc StarkCurve) Add(x1, y1, x2, y2 *big.Int) (x, y *big.Int) {
	x, y, z := sc.add(x1, y1, big.NewInt(1), x2, y2, big.NewInt(1))
	return DivMod(x, z, sc.P), DivMod(y, z, sc.P)
}

func (sc StarkCurve) add(x1, y1, z1, x2, y2, z2 *big.Int) (x, y, z *big.Int) {
	// As elliptic curves form a group, there is an additive identity that is the equivalent of 0
	// If 𝑃=0 or 𝑄=0, then 𝑃+𝑄=𝑄 or 𝑃+𝑄=𝑃, respectively
	// NOTICE: the EC multiplication algorithm is using using `StarkCurve.rewriteScalar` trick
	//   to avoid this condition and provide constant-time execution.

	if (len(x1.Bits()) == 0 && len(y1.Bits()) == 0) || len(z1.Bits()) == 0 {
		return x2, y2, z2
	}
	if (len(x2.Bits()) == 0 && len(y2.Bits()) == 0) || len(z2.Bits()) == 0 {
		return x1, y1, z1
	}

	// 1. t0 = X1 * X2
	t0 := new(big.Int).Mul(x1, x2)
	t0.Mod(t0, sc.P)

	// 2. t1 = Y1 * Y2
	t1 := new(big.Int).Mul(y1, y2)
	t1.Mod(t1, sc.P)

	// 3. t3 = X2 + Y2
	t3 := new(big.Int).Add(x2, y2)
	t3.Mod(t3, sc.P)

	// 4. t4 = X1 + Y1
	t4 := new(big.Int).Add(x1, y1)
	t4.Mod(t4, sc.P)

	// 5. t3 = t3 * t4
	t3.Mul(t3, t4)
	t3.Mod(t3, sc.P)

	// 6. t4 = t0 + t1
	t4.Add(t0, t1)
	t4.Mod(t4, sc.P)

	// 7. t3 = t3 - t4
	t3.Sub(t3, t4)
	t3.Mod(t3, sc.P)

	// 8. t4 = X2 * Z1
	t4.Mul(x2, z1)
	t4.Mod(t4, sc.P)

	// 9. t4 = t4 + X1
	t4.Add(t4, x1)
	t4.Mod(t4, sc.P)

	// 10. t5 = Y2 * Z1
	t5 := new(big.Int).Mul(y2, z1)
	t5.Mod(t5, sc.P)

	// 11. t5 = t5 + Y1
	t5.Add(t5, y1)
	t5.Mod(t5, sc.P)

	// 12. Z3 = a * t4
	z3 := new(big.Int).Mul(sc.Alpha, t4)
	z3.Mod(z3, sc.P)

	// 13. X3 = b3 * Z1
	x3 := new(big.Int).Mul(Curve.B3, z1)
	x3.Mod(x3, sc.P)

	// 14. Z3 = X3 + Z3
	z3.Add(x3, z3)
	z3.Mod(z3, sc.P)

	// 15. X3 = t1 - Z3
	x3.Sub(t1, z3)
	x3.Mod(x3, sc.P)

	// 16. Z3 = t1 + Z3
	z3.Add(t1, z3)
	z3.Mod(z3, sc.P)

	// 17. Y3 = X3 * Z3
	y3 := new(big.Int).Mul(x3, z3)
	y3.Mod(y3, sc.P)

	// 18. t1 = t0 + t0
	t1.Add(t0, t0)
	t1.Mod(t1, sc.P)

	// 19. t1 = t1 + t0
	t1.Add(t1, t0)
	t1.Mod(t1, sc.P)

	// 20. t2 = a * Z1
	t2 := new(big.Int).Mul(sc.Alpha, z1)
	t2.Mod(t2, sc.P)

	// 21. t4 = b3 * t4
	t4.Mul(Curve.B3, t4)
	t4.Mod(t4, sc.P)

	// 22. t1 = t1 + t2
	t1.Add(t1, t2)
	t1.Mod(t1, sc.P)

	// 23. t2 = t0 - t2
	t2.Sub(t0, t2)
	t2.Mod(t2, sc.P)

	// 24. t2 = a * t2
	t2.Mul(sc.Alpha, t2)
	t2.Mod(t2, sc.P)

	// 25. t4 = t4 + t2
	t4.Add(t4, t2)
	t4.Mod(t4, sc.P)

	// 26. t0 = t1 * t4
	t0.Mul(t1, t4)
	t0.Mod(t0, sc.P)

	// 27. Y3 = Y3 + t0
	y3.Add(y3, t0)
	y3.Mod(y3, sc.P)

	// 28. t0 = t5 * t4
	t0.Mul(t5, t4)
	t0.Mod(t0, sc.P)

	// 29. X3 = t3 * X3
	x3.Mul(t3, x3)
	x3.Mod(x3, sc.P)

	// 30. X3 = X3 - t0
	x3.Sub(x3, t0)
	x3.Mod(x3, sc.P)

	// 31. t0 = t3 * t1
	t0.Mul(t3, t1)
	t0.Mod(t0, sc.P)

	// 32. Z3 = t5 * Z3
	z3.Mul(t5, z3)
	z3.Mod(z3, sc.P)

	// 33. Z3 = Z3 + t0
	z3.Add(z3, t0)
	z3.Mod(z3, sc.P)

	return x3, y3, z3
}

// Doubles a point on an elliptic curve with the equation y^2 = x^3 + alpha*x + beta mod p.
// Assumes affine form (x, y) is spread (x1 *big.Int, y1 *big.Int)
//
// This implements Algorithm 3 of 2015 Renes–Costello–Batina "Complete addition formulas for prime order elliptic curves"
// (ref: https://eprint.iacr.org/2015/1060.pdf)
func (sc StarkCurve) Double(x1, y1 *big.Int) (x_out, y_out *big.Int) {
	x, y, z := sc.double(x1, y1, big.NewInt(1))
	return DivMod(x, z, sc.P), DivMod(y, z, sc.P)
}

func (sc StarkCurve) double(x1, y1, z1 *big.Int) (x_out, y_out, z_out *big.Int) {
	// 1. t0 = X * X
	t0 := new(big.Int).Mul(x1, x1)
	t0.Mod(t0, sc.P)

	// 2. t1 = Y * Y
	t1 := new(big.Int).Mul(y1, y1)
	t1.Mod(t1, sc.P)

	// 3. t2 = Z * Z
	t2 := new(big.Int).Mul(z1, z1)
	t2.Mod(t2, sc.P)

	// 4. t3 = X * Y
	t3 := new(big.Int).Mul(x1, y1)
	t3.Mod(t3, sc.P)

	// 5. t3 = t3 + t3
	t3.Add(t3, t3)
	t3.Mod(t3, sc.P)

	// 6. Z3 = X * Z
	z3 := new(big.Int).Mul(x1, z1)
	z3.Mod(z3, sc.P)

	// 7. Z3 = Z3 + Z3
	z3.Add(z3, z3)
	z3.Mod(z3, sc.P)

	// 8. X3 = a * Z3
	x3 := new(big.Int).Mul(sc.Alpha, z3)
	x3.Mod(x3, sc.P)

	// 9. Y3 = b3 * t2
	y3 := new(big.Int).Mul(Curve.B3, t2)
	y3.Mod(y3, sc.P)

	// 10. Y3 = X3 + Y3
	y3.Add(x3, y3)
	y3.Mod(y3, sc.P)

	// 11. X3 = t1 - Y3
	x3.Sub(t1, y3)
	x3.Mod(x3, sc.P)

	// 12. Y3 = t1 + Y3
	y3.Add(t1, y3)
	y3.Mod(y3, sc.P)

	// 13. Y3 = X3 * Y3
	y3.Mul(x3, y3)
	y3.Mod(y3, sc.P)

	// 14. X3 = t3 * X3
	x3.Mul(t3, x3)
	x3.Mod(x3, sc.P)

	// 15. Z3 = b3 * Z3
	z3.Mul(Curve.B3, z3)
	z3.Mod(z3, sc.P)

	// 16. t2 = a * t2
	t2.Mul(sc.Alpha, t2)
	t2.Mod(t2, sc.P)

	// 17. t3 = t0 - t2
	t3.Sub(t0, t2)
	t3.Mod(t3, sc.P)

	// 18. t3 = a * t3
	t3.Mul(sc.Alpha, t3)
	t3.Mod(t3, sc.P)

	// 19. t3 = t3 + z3
	t3.Add(t3, z3)
	t3.Mod(t3, sc.P)

	// 20. z3 = t0 + t0
	z3.Add(t0, t0)
	z3.Mod(z3, sc.P)

	// 21. t0 = Z3 + t0
	t0.Add(z3, t0)
	t0.Mod(t0, sc.P)

	// 22. t0 = t0 + t2
	t0.Add(t0, t2)
	t0.Mod(t0, sc.P)

	// 23. t0 = t0 * t3
	t0.Mul(t0, t3)
	t0.Mod(t0, sc.P)

	// 24. Y3 = Y3 + t0
	y3.Add(y3, t0)
	y3.Mod(y3, sc.P)

	// 25. t2 = Y * Z
	t2.Mul(y1, z1)
	t2.Mod(t2, sc.P)

	// 26. t2 = t2 + t2
	t2.Add(t2, t2)
	t2.Mod(t2, sc.P)

	// 27. t0 = t2 * t3
	t0.Mul(t2, t3)
	t0.Mod(t0, sc.P)

	// 28. X3 = X3 - t0
	x3.Sub(x3, t0)
	x3.Mod(x3, sc.P)

	// 29. Z3 = t2 * t1
	z3.Mul(t2, t1)
	z3.Mod(z3, sc.P)

	// 30. Z3 = Z3 + Z3
	z3.Add(z3, z3)
	z3.Mod(z3, sc.P)

	// 31. Z3 = Z3 + Z3
	z3.Add(z3, z3)
	z3.Mod(z3, sc.P)

	return x3, y3, z3
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
	// Infinity is considered part of the curve
	if len(x.Bits()) == 0 && len(y.Bits()) == 0 {
		return true
	}

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
	var _ecMult = func(m, x1, y1, z1 *big.Int) (x, y, z *big.Int) {
		// Two-index table initialization, Q[0] <- P
		q := [2]struct {
			x *big.Int
			y *big.Int
			z *big.Int
		}{
			{
				x: x1,
				y: y1,
				z: z1,
			},
			{
				x: nil,
				y: nil,
				z: nil,
			},
		}

		// Run the algorithm, expects the most-significant bit is 1
		for i := sc.N.BitLen() - 2; i >= 0; i-- {
			q[0].x, q[0].y, q[0].z = sc.double(q[0].x, q[0].y, q[0].z)          // Q[0] <- 2Q[0]
			q[1].x, q[1].y, q[1].z = sc.add(q[0].x, q[0].y, q[0].z, x1, y1, z1) // Q[1] <- Q[0] + P
			b := m.Bit(i)                                                       // b    <- bit at position i
			q[0].x, q[0].y, q[0].z = q[b].x, q[b].y, q[b].z                     // Q[0] <- Q[b]
		}

		return q[0].x, q[0].y, q[0].z
	}

	xOut, yOut, zOut := _ecMult(sc.rewriteScalar(m), x1, y1, big.NewInt(1))
	return DivMod(xOut, zOut, sc.P), DivMod(yOut, zOut, sc.P)
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
	mod := new(big.Int).Mod(size, sc.N)                       // 2ˆn mod q
	diff := new(big.Int).Sub(k, mod)                          // (k - 2ˆn mod q)
	return new(big.Int).Add(size, diff)                       // 2ˆn + (k - 2ˆn mod q)
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
