package caigo

import (
	"fmt"
	"time"
	"math"
	"math/big"
	"testing"

	"gonum.org/v1/gonum/stat"
)

func BenchmarkPedersenHash(b *testing.B) {
	suite := [][]*big.Int{
		{HexToBN("0x12773"), HexToBN("0x872362")},
		{HexToBN("0x1277312773"), HexToBN("0x872362872362")},
		{HexToBN("0x1277312773"), HexToBN("0xCD2a3d9F938E13CD947Ec05AbC7FE734Df8DD826")},
		{HexToBN("0xbBbBBBBbbBBBbbbBbbBbbbbBBbBbbbbBbBbbBBbB"), HexToBN("0x872362872362")},
		{HexToBN("0xCD2a3d9F938E13CD947Ec05AbC7FE734Df8DD826"), HexToBN("0xbBbBBBBbbBBBbbbBbbBbbbbBBbBbbbbBbBbbBBbB")},
		{HexToBN("0x7f15c38ea577a26f4f553282fcfe4f1feeb8ecfaad8f221ae41abf8224cbddd"), HexToBN("0x13d41f388b8ea4db56c5aa6562f13359fab192b3db57651af916790f9debee9")},
		{HexToBN("0x7f15c38ea577a26f4f553282fcfe4f1feeb8ecfaad8f221ae41abf8224cbddd"), HexToBN("0x7f15c38ea577a26f4f553282fcfe4f1feeb8ecfaad8f221ae41abf8224cbdde")},
	}

	for _, test := range suite {
		b.Run(fmt.Sprintf("input_size_%d_%d", test[0].BitLen(), test[1].BitLen()), func(b *testing.B) {
			Curve.PedersenHash(test)
		})
	}
}

func TestPedersenHash(t *testing.T) {
	testPedersen := []struct {
		elements []*big.Int
		expected *big.Int
	}{
		{
			elements: []*big.Int{HexToBN("0x12773"), HexToBN("0x872362")},
			expected: HexToBN("0x5ed2703dfdb505c587700ce2ebfcab5b3515cd7e6114817e6026ec9d4b364ca"),
		},
		{
			elements: []*big.Int{HexToBN("0x13d41f388b8ea4db56c5aa6562f13359fab192b3db57651af916790f9debee9"), HexToBN("0x537461726b4e6574204d61696c")},
			expected: HexToBN("0x180c0a3d13c1adfaa5cbc251f4fc93cc0e26cec30ca4c247305a7ce50ac807c"),
		},
		{
			elements: []*big.Int{big.NewInt(100), big.NewInt(1000)},
			expected: HexToBN("0x45a62091df6da02dce4250cb67597444d1f465319908486b836f48d0f8bf6e7"),
		},
	}

	for _, tt := range testPedersen {
		hash, err := Curve.PedersenHash(tt.elements)
		if err != nil {
			t.Errorf("Hashing err: %v\n", err)
		}
		if hash.Cmp(tt.expected) != 0 {
			t.Errorf("incorrect hash: got %v expected %v\n", hash, tt.expected)
		}
	}
}

func TestDivMod(t *testing.T) {
	testDivmod := []struct {
		x        *big.Int
		y        *big.Int
		expected *big.Int
	}{
		{
			x:        StrToBig("311379432064974854430469844112069886938521247361583891764940938105250923060"),
			y:        StrToBig("621253665351494585790174448601059271924288186997865022894315848222045687999"),
			expected: StrToBig("2577265149861519081806762825827825639379641276854712526969977081060187505740"),
		},
		{
			x:        big.NewInt(1),
			y:        big.NewInt(2),
			expected: HexToBN("0x0400000000000008800000000000000000000000000000000000000000000001"),
		},
	}

	for _, tt := range testDivmod {
		divR := DivMod(tt.x, tt.y, Curve.P)

		if divR.Cmp(tt.expected) != 0 {
			t.Errorf("DivMod Res %v does not == expected %v\n", divR, tt.expected)
		}
	}
}

func TestAdd(t *testing.T) {
	testAdd := []struct {
		x         *big.Int
		y         *big.Int
		expectedX *big.Int
		expectedY *big.Int
	}{
		{
			x:         StrToBig("1468732614996758835380505372879805860898778283940581072611506469031548393285"),
			y:         StrToBig("1402551897475685522592936265087340527872184619899218186422141407423956771926"),
			expectedX: StrToBig("2573054162739002771275146649287762003525422629677678278801887452213127777391"),
			expectedY: StrToBig("3086444303034188041185211625370405120551769541291810669307042006593736192813"),
		},
		{
			x:         big.NewInt(1),
			y:         big.NewInt(2),
			expectedX: StrToBig("225199957243206662471193729647752088571005624230831233470296838210993906468"),
			expectedY: StrToBig("190092378222341939862849656213289777723812734888226565973306202593691957981"),
		},
	}

	for _, tt := range testAdd {
		resX, resY := Curve.Add(Curve.Gx, Curve.Gy, tt.x, tt.y)
		if resX.Cmp(tt.expectedX) != 0 {
			t.Errorf("ResX %v does not == expected %v\n", resX, tt.expectedX)

		}
		if resY.Cmp(tt.expectedY) != 0 {
			t.Errorf("ResY %v does not == expected %v\n", resY, tt.expectedY)
		}
	}
}

func TestMultAir(t *testing.T) {
	testMult := []struct {
		r         *big.Int
		x         *big.Int
		y         *big.Int
		expectedX *big.Int
		expectedY *big.Int
	}{
		{
			r:         StrToBig("2458502865976494910213617956670505342647705497324144349552978333078363662855"),
			x:         StrToBig("1468732614996758835380505372879805860898778283940581072611506469031548393285"),
			y:         StrToBig("1402551897475685522592936265087340527872184619899218186422141407423956771926"),
			expectedX: StrToBig("182543067952221301675635959482860590467161609552169396182763685292434699999"),
			expectedY: StrToBig("3154881600662997558972388646773898448430820936643060392452233533274798056266"),
		},
	}

	for _, tt := range testMult {
		x, y, err := Curve.MimicEcMultAir(tt.r, tt.x, tt.y, Curve.Gx, Curve.Gy)
		if err != nil {
			t.Errorf("MultAirERR %v\n", err)
		}

		if x.Cmp(tt.expectedX) != 0 {
			t.Errorf("ResX %v does not == expected %v\n", x, tt.expectedX)

		}
		if y.Cmp(tt.expectedY) != 0 {
			t.Errorf("ResY %v does not == expected %v\n", y, tt.expectedY)
		}
	}
}


func TestEcMult(t *testing.T) {
	testMult := []struct {
		k         *big.Int
		expectedX *big.Int
		expectedY *big.Int
	}{
		{
			k:         StrToBig("1"),
			expectedX: StrToBig("874739451078007766457464989774322083649278607533249481151382481072868806602"),
			expectedY: StrToBig("152666792071518830868575557812948353041420400780739481342941381225525861407"),
		},
		{
			k:         StrToBig("2"),
			expectedX: StrToBig("3324833730090626974525872402899302150520188025637965566623476530814354734325"),
			expectedY: StrToBig("3147007486456030910661996439995670279305852583596209647900952752170983517249"),
		},
		{
			k:         StrToBig("3"),
			expectedX: StrToBig("1839793652349538280924927302501143912227271479439798783640887258675143576352"),
			expectedY: StrToBig("3564972295958783757568195431080951091358810058262272733141798511604612925062"),
		},
	}

	for _, tt := range testMult {
		x, y, err := Curve.PrivateToPoint(tt.k)

		if err != nil {
			t.Errorf("EcMult %v\n", err)
		}

		if x.Cmp(tt.expectedX) != 0 {
			t.Errorf("ResX %v does not == expected %v\n", x, tt.expectedX)

		}
		if y.Cmp(tt.expectedY) != 0 {
			t.Errorf("ResY %v does not == expected %v\n", y, tt.expectedY)
		}
	}
}

func BenchmarkPrivateToPoint(b *testing.B) {
	var _genNBits = func(n int) (k *big.Int) {
		k = big.NewInt(1)
		for i := 0; i < n; i++ {
			k = k.Lsh(k, 1).Add(k, big.NewInt(1))
		}
		return
	}

	xs := []float64{}
	for i := 1; i < Curve.N.BitLen() - 1; i++ {
		k := _genNBits(i)
		b.Run(fmt.Sprintf("input_size_%d", k.BitLen()), func(b *testing.B) {
			start := time.Now()
			Curve.PrivateToPoint(k)
			elapsed := time.Since(start).Nanoseconds()
			xs = append(xs, float64(elapsed))
		})
	}

	// computes the weighted mean of the dataset.
	// we don't have any weights (ie: all weights are 1)
	// so we just pass a nil slice.
	mean := stat.Mean(xs, nil)
	variance := stat.Variance(xs, nil)
	stddev := math.Sqrt(variance)
	fmt.Printf("mean=     %v\n", mean)
	fmt.Printf("variance= %v\n", variance)
	fmt.Printf("std-dev=  %v\n", stddev)
}
