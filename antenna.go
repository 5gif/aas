package aas

import (
	"math"
	"math/cmplx"

	log "github.com/sirupsen/logrus"
	"github.com/wiless/vlib"
)

type Antenna struct {
	AntennaConfig   []int     `json:"AntennaConfig"`
	ElectricalTilt  []float64 `json:"ElectricalTilt"`
	Escan           []float64 `json:"Escan"`
	GainDb          float64   `json:"GainDb"`
	HBeamWidth      float64   `json:"HBeamWidth"`
	Omni            bool      `json:"Omni"`
	SLAV            float64   `json:"SLAV"`
	VBeamWidth      float64   `json:"VBeamWidth"`
	EspacingHfactor float64   `json:"EspacingHfactor"`
	EspacingVfactor float64   `json:"EspacingVfactor"`
	MechanicalTilt  float64   `json:"MechanicalTilt"`
	Polarization    []float64 `json:"Polarization"`
	PanelAz         []float64 `json:"PanelAz"`
	PanelEl         []float64 `json:"PanelEl"`
}

func (ant *Antenna) GaindB(theta, phi float64) (aag map[int]vlib.MatrixF, bestBeamID int, Az, El float64) {
	theta = wrap180To180(theta)
	phi = wrap0To180(phi)
	var ag float64
	Az, El, ag = ant.ElementGainDb(theta, phi)
	hspace := ant.EspacingHfactor
	vspace := ant.EspacingVfactor
	var sum = complex(0.0, 0.0)

	dtilt := ant.ElectricalTilt // degree
	descan := ant.Escan         //degree

	nv := ant.AntennaConfig[0] / ant.AntennaConfig[5]
	nh := ant.AntennaConfig[1] / ant.AntennaConfig[6]

	maxgain := -1000.0
	bestBeamID = 0
	nbeams := len(ant.Escan) * len(ant.ElectricalTilt)
	aag = make(map[int]vlib.MatrixF, nbeams)
	var c = complex(math.Sqrt(1/float64(nv*nh)), 0)
	for i := 0; i < len(dtilt); i++ { //  dtilt is a vector of Zenith Angles of the Beam Set
		for j := 0; j < len(descan); j++ { // descan is a vector of Azimuth Angles of the Beam Set
			beamid := j + len(descan)*i
			sum = 0.0
			phiP := -math.Cos(dtilt[i]*math.Pi/180) + math.Cos((phi-ant.MechanicalTilt+90)*math.Pi/180)
			phiR := -math.Sin(dtilt[i]*math.Pi/180)*math.Sin(descan[j]*math.Pi/180) + math.Sin((phi-ant.MechanicalTilt+90)*math.Pi/180)*math.Sin(theta*math.Pi/180)
			for m := 1; m <= nv; m++ {
				for n := 1; n <= nh; n++ {
					w := cmplx.Exp(complex(0, 2*math.Pi*(float64(m-1)*vspace*phiP)))
					v := cmplx.Exp(complex(0, 2*math.Pi*(float64(n-1)*hspace*phiR)))
					sum = sum + c*w*v
				}
			}
			txRUGains := vlib.NewMatrixF(ant.AntennaConfig[5], ant.AntennaConfig[6])
			for k := 0; k < ant.AntennaConfig[5]; k++ {
				for l := 0; l < ant.AntennaConfig[6]; l++ {
					txRUGains[k][l] = ag + (10 * math.Log10(math.Pow(cmplx.Abs(sum), 2))) // Composite Beam Gain + Antenna Element Gain
					temp := txRUGains[k][l]
					if maxgain < temp {
						maxgain = temp
						bestBeamID = beamid
					}
				}
			}
			aag[beamid] = txRUGains
		}
	}
	return aag, bestBeamID, Az, El

}

func (ant *Antenna) ElementGainDb(theta, phi float64) (az, el, Ag float64) {
	phi = wrap0To180(phi)
	theta = wrap180To180(theta)
	MaxGaindBi := ant.GainDb //    0 for ue and 8 for bs
	theta3dB := ant.GainDb   // degree
	phi3dB := ant.VBeamWidth
	SLAmax := ant.SLAV
	Am := SLAmax
	Ah := -math.Min(12.0*math.Pow(theta/theta3dB, 2.0), Am)
	MechTiltGCS := ant.MechanicalTilt // Pointing to Horizon..axis..
	Av := -math.Min(12.0*math.Pow((phi-MechTiltGCS)/phi3dB, 2.0), SLAmax)
	result := -math.Min(-math.Floor(Av+Ah), Am)
	//result = Ah
	az = Ah
	el = Av
	Ag = result + MaxGaindBi
	return az, el, Ag
}

func (ant *Antenna) GetPorts() int {
	//polarization can be 1 or 2
	//error handling to be done

	polar := ant.AntennaConfig[2]
	pv := ant.AntennaConfig[3]
	ph := ant.AntennaConfig[4]
	if pv > 1 || ph > 1 {
		log.Panic("Unsupported Panels >1")
	}

	p := ant.AntennaConfig[5] * ant.AntennaConfig[6] * polar * pv * ph
	return p
}

func (ant *Antenna) FindTxruLocation() []vlib.MatrixF {
	//Done only for 1 panel
	if ant.AntennaConfig[2]*ant.AntennaConfig[3] > 1 {
		log.Panic("Number of Panels more than 1")
	}

	var A, B, C, D vlib.Location3D
	var Centre [][]vlib.Location3D
	p := ant.GetPorts()

	N := float64(ant.AntennaConfig[0])
	M := float64(ant.AntennaConfig[1])
	Np := float64(ant.AntennaConfig[5])
	Mp := float64(ant.AntennaConfig[6])
	Centre = make([][]vlib.Location3D, int(Np))
	dh := ant.EspacingHfactor
	dv := ant.EspacingVfactor

	A.X = 0.0
	A.Y = 0.0
	A.Z = 0.0
	B.X = 0.0
	B.Y = (N * dv) / Np
	B.Z = 0.0
	C.X = (M * dh) / Mp
	C.Y = (N * dv) / Np
	C.Z = 0.0
	D.X = (M * dh) / Mp
	D.Y = 0.0
	D.Z = 0.0

	var ref vlib.Location3D
	ref.X = D.X - A.X
	ref.Y = B.Y - A.Y
	ref.Z = B.Z - A.Z

	for i := 0; i < int(Np); i++ {
		Centre[i] = make([]vlib.Location3D, int(Mp))

		Centre[i][0].Y = ref.Y + dv*Np*float64(i)

		for j := 0; j < int(Mp); j++ {
			Centre[i][j].X = ref.X + dh*Mp*float64(j)
			Centre[i][j].Y = Centre[i][0].Y
			Centre[i][j].Z = ref.Z

		}
	}

	var Dx = make([]vlib.MatrixF, p)
	ind := 0

	for ii := 0; ii < ant.AntennaConfig[2]; ii++ {
		for i := 0; i < int(Np); i++ {
			for j := 0; j < int(Mp); j++ {
				drx := vlib.NewMatrixF(3, 1)
				drx[0][0] = Centre[i][j].X
				drx[1][0] = Centre[i][j].Y
				drx[2][0] = Centre[i][j].Z
				Dx[ind] = drx
				ind = ind + 1

			}

		}

	}

	return Dx

}

// Wrap0To180 wraps the input angle to 0 to 180
func wrap0To180(degree float64) float64 {
	if degree >= 0 && degree <= 180 {
		return degree
	}
	if degree < 0 {
		degree = -degree
	}
	if degree >= 360 {
		degree = math.Mod(degree, 360)
	}
	if degree > 180 {

		degree = 360 - degree
	}
	return degree
}

// Wrap180To180 wraps the input angle to -180 to 180
func wrap180To180(degree float64) float64 {
	if degree >= -180 && degree <= 180 {
		return degree
	}
	if degree > 180 {
		rem := math.Mod(degree, 180.0)
		degree = -180 + rem

	} else if degree < -180 {
		rem := math.Mod(degree, 180.0)
		//	fmt.Println("Remainder for ", degree, rem)
		degree = 180 + rem
	}
	return degree
}
