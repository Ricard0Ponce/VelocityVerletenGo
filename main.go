package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

type Celda struct {
	nac  int
	iats []int
}
type Atomo struct {
	P    []float64 // Posiciones
	V    []float64 // Velocidades
	A    []float64 // Aceleraciones
	NV   int       // Número de vecinos de cada átomo
	VC   []int     // Vector de vecinos cercanos
	NCel int       // Número de celda
	CCel [3]int    // Coordenadas de la celda
}

type Sim struct {
	temp, dens, L, d, r, rc, v0, drc, rcd, dl           float64 // nrc
	naa, nd, na, nc, ncL, nct, ncvc, ncp, nci, icc, nrc int
	cuad, metodo                                        string
	dt                                                  float64
	atoms                                               []*Atomo
	dis                                                 []float64
	f                                                   float64
	et, ec, ei                                          float64
	celdas                                              []*Celda
	nctt                                                [3]int
}

func (s *Sim) Datos() {
	fmt.Println()
	fmt.Println("Sim::Datos")

	s.temp = 1.0
	s.dens = 0.65

	s.naa = 15                                          // Numero de atomos por lado // 8
	s.nd = 3                                            // Numero de dimensiones
	s.na = int(math.Pow(float64(s.naa), float64(s.nd))) // Numero de atomos
	s.nc = 1000                                         // Numero de configuraciones
	s.ncL = int(s.L)
	s.nct = 0
	s.cuad = "cuad"
	//s.metodo = "orig"
	//s.metodo = "vc"
	s.metodo = "cel"
	s.dt = 0.001 // Delta de tiempo
	s.d = 1.0
	s.L = math.Pow(float64(s.na)/s.dens, 1.0/float64(s.nd))
	s.r = s.d / 2.0
	s.rc = 3.0
	s.v0 = 1.2 * math.Sqrt(s.temp*float64(s.nd)/2.0)
	s.drc = 0.5
	s.rcd = s.rc + s.drc
	s.ncvc = 10
	s.nrc = int(s.rc)
	s.ncp = 10 // Numero de configuraciones parciales
	s.nci = s.nc / s.ncp
	s.dis = make([]float64, s.nd)
}

func (s *Sim) EscribirDatos() {
	fmt.Println()
	fmt.Println("Sim::EscribirDatos")

	fmt.Println()
	fmt.Println("Datos iniciales")
	fmt.Printf("na,nd,nc %d %d %d\n", s.na, s.nd, s.nc)
	fmt.Printf("temp,dens %.2f %.2f\n", s.temp, s.dens)
	fmt.Printf("dt,L,d,r %.3f %.3f %.3f %.3f\n", s.dt, s.L, s.d, s.r)
	fmt.Printf("cuad, met %s %s\n", s.cuad, s.metodo)
	fmt.Printf("v0, rc %.3f %.3f\n", s.v0, s.rc)
	fmt.Printf("ncp,nci %d %d\n", s.ncp, s.nci)
}

func (s *Sim) IniciarAtomos() {
	fmt.Println()
	fmt.Println("Sim::IniciarAtomos")

	pv := make([]float64, 3)

	for ia := 0; ia < s.na; ia++ {
		att := &Atomo{}
		s.atoms = append(s.atoms, att)
		//s.atoms[ia].A = append(s.atoms[ia].A, 0.0)
	}
	an := len(s.atoms)
	fmt.Println("an", an)

	for ia := 0; ia < s.na; ia++ {
		for id := 0; id < s.nd; id++ {
			pv[id] = s.L * (rand.Float64())
			s.atoms[ia].P = append(s.atoms[ia].P, pv[id])
			s.atoms[ia].A = append(s.atoms[ia].A, 0.0)
		}
	}

	for ia := 0; ia < s.na; ia++ {
		for id := 0; id < s.nd; id++ {
			r := rand.Float64()
			r = (r * 2.0) - 1
			pv[id] = s.v0 * r
			s.atoms[ia].V = append(s.atoms[ia].V, pv[id])
		}
	}
}

func (s *Sim) Aceleraciones() {
	for ia := 0; ia < s.na; ia++ {
		for id := 0; id < s.nd; id++ {
			s.atoms[ia].A[id] = 0.0
		}
	}

	var u float64
	for ia := 0; ia < s.na-1; ia++ {
		for ja := ia + 1; ja < s.na; ja++ {
			r := s.Dist(ia, ja)
			u += s.LJ(r)

			for id := 0; id < s.nd; id++ {
				s.atoms[ia].A[id] += s.f * s.dis[id]
				s.atoms[ja].A[id] -= s.f * s.dis[id]
			}
		}
	}
}

func (s *Sim) LJ(r float64) float64 {
	r2 := 1.0 / r
	r6 := r2 * r2 * r2
	r12 := r6 * r6

	u := 4.0 * (r12 - r6)
	s.f = 24.0 * (2.0*r12 - r6) * r2

	return u
}

func (s *Sim) Dist(i, j int) float64 {
	var r float64

	for id := 0; id < s.nd; id++ {
		s.dis[id] = s.atoms[i].P[id] - s.atoms[j].P[id]
		if math.Abs(s.dis[id]) > s.L/2 {
			s.dis[id] -= math.Copysign(s.L, s.dis[id])
		}
		r += s.dis[id] * s.dis[id]
	}

	return r
}

func (s *Sim) Prop() {
	s.et = 0
	s.ec = 0
	s.ei = 0

	for ia := 0; ia < s.na; ia++ {
		for id := 0; id < s.nd; id++ {
			vi := s.atoms[ia].V[id]
			s.ec += (vi * vi) * 0.5
		}
	}
	s.ec /= float64(s.na)

	for ia := 0; ia < s.na-1; ia++ {
		for ib := ia + 1; ib < s.na; ib++ {
			r := s.Dist(ia, ib)
			s.ei += s.LJ(r)
		}
	}
	s.ei /= float64(s.na)
	s.et = s.ec + s.ei
	s.temp = (s.ec * 2.0) / float64(s.nd)
}

func (s *Sim) Cuadrada() {
	fmt.Println()
	fmt.Println("Sim::Cuadrada")

	nl := int(s.L)
	dl := s.L / float64(nl)
	var x, y, z float64
	ia := 0

	for iz := 0; iz < s.naa; iz++ {
		z = float64(iz) * dl

		for iy := 0; iy < s.naa; iy++ {
			y = float64(iy) * dl
			for ix := 0; ix < s.naa; ix++ {
				x = float64(ix) * dl
				s.atoms[ia].P = []float64{x, y, z}
				ia++
			}
		}

		if s.nd == 2 {
			break
		}
	}
}

func (s *Sim) calcularVC() {
	for ia := 0; ia < s.na; ia++ {
		s.atoms[ia].VC = nil
		s.atoms[ia].NV = 0
	}

	for ia := 0; ia < s.na-1; ia++ {
		for ja := ia + 1; ja < s.na; ja++ {
			r := s.Dist(ia, ja) // Se calculo la distancia
			if r < s.rcd {
				s.atoms[ia].VC = append(s.atoms[ia].VC, ja)
				s.atoms[ja].VC = append(s.atoms[ja].VC, ia)
			}
		}
		s.atoms[ia].NV = len(s.atoms[ia].VC)
	}
	s.atoms[s.na-1].NV = len(s.atoms[s.na-1].VC)
}

func (s *Sim) AceleracionesVC() {
	for ia := 0; ia < s.na; ia++ {
		for id := 0; id < s.nd; id++ {
			s.atoms[ia].A[id] = 0.0
		}
	}

	var u float64
	for ia := 0; ia < s.na; ia++ {
		nv := s.atoms[ia].NV
		for ja := 0; ja < nv; ja++ {
			ja1 := s.atoms[ia].VC[ja]
			r := s.Dist(ia, ja1) // Se calculo la distancia
			if r < s.rc {
				u += s.LJ(r)
				for id := 0; id < s.nd; id++ {
					s.atoms[ia].A[id] += s.f * s.dis[id] * 0.5
					// s.atoms[ja1].A[id] -= s.f * s.dis[id] / r
				}
			}
		}
	}
}

func (s *Sim) CalcularCel() {
	var icel, ia, id, acel, api int
	var ap float64

	// Ensure celdas is initialized
	//if s.celdas == nil || len(s.celdas) == 0 {
	//	s.IniciarCel()
	//}

	for icel = 0; icel < s.nct; icel++ {
		s.celdas[icel].nac = 0
		s.celdas[icel].iats = nil
	}

	for ia = 0; ia < s.na; ia++ {
		acel = 0

		for id = 0; id < s.nd; id++ {
			ap = s.atoms[ia].P[id]
			api = int(ap / s.dl)
			if api < 0 {
				api += s.ncL
			}
			if api >= s.ncL {
				api -= s.ncL
			}
			//if api < 0 || api >= s.ncL {
			//	panic(fmt.Sprintf("api out of range: %d", api))
			//}
			s.atoms[ia].CCel[id] = api
			acel += api * s.nctt[id]
		}

		//if acel < 0 || acel >= s.nct {
		//	panic(fmt.Sprintf("acel out of range: %d", acel))
		//}

		s.atoms[ia].NCel = acel
		s.celdas[acel].nac++
		s.celdas[acel].iats = append(s.celdas[acel].iats, ia)
	}
}

func (s *Sim) IniciarCel() {
	fmt.Println()
	fmt.Println("Sim::InicCel")

	s.ncL = int(s.L)
	s.dl = s.L / float64(s.ncL)
	s.nct = s.ncL * s.ncL
	if s.nd == 3 {
		s.nct *= s.ncL
	}

	s.celdas = make([]*Celda, s.nct)
	for i := 0; i < s.nct; i++ {
		s.celdas[i] = &Celda{}
		s.celdas[i].nac = 0
	}

	s.nctt[0] = 1
	s.nctt[1] = s.ncL
	s.nctt[2] = s.ncL * s.ncL
}

func (s *Sim) AceleracionesCel() {
	for ia := 0; ia < s.na; ia++ {
		for id := 0; id < s.nd; id++ {
			s.atoms[ia].A[id] = 0.0
		}
	}

	var u float64 = 0
	for ia := 0; ia < s.na; ia++ {
		for iz := -s.nrc; iz <= s.nrc; iz++ {
			z := float64(iz) + float64(s.atoms[ia].CCel[2])
			if z < 0 {
				z += float64(s.ncL)
			}
			if z > float64(s.ncL-1) {
				z -= float64(s.ncL)
			}

			for iy := -s.nrc; iy <= s.nrc; iy++ {
				y := float64(iy) + float64(s.atoms[ia].CCel[1])
				if y < 0 {
					y += float64(s.ncL)
				}
				if y > float64(s.ncL-1) {
					y -= float64(s.ncL)
				}

				for ix := -s.nrc; ix <= s.nrc; ix++ {
					x := float64(ix) + float64(s.atoms[ia].CCel[0])
					if x < 0 {
						x += float64(s.ncL)
					}
					if x > float64(s.ncL-1) {
						x -= float64(s.ncL)
					}
					ncel := int(x)*s.nctt[0] + int(y)*s.nctt[1] + int(z)*s.nctt[2]
					nac := s.celdas[ncel].nac

					if nac > 0 {
						for ja := 0; ja < nac; ja++ {
							jal := s.celdas[ncel].iats[ja]
							if ia > jal {
								r := s.Dist(ia, jal)
								if r < s.rc*s.rc {
									u += s.LJ(r)
									for id := 0; id < s.nd; id++ {
										s.atoms[ia].A[id] += s.f * s.dis[id]
										s.atoms[jal].A[id] -= s.f * s.dis[id]
									}
								}
							}
						}
					}
				}
			}
		}
	}
	u /= float64(s.na)
}

func (s *Sim) Simulacion() {
	fmt.Println()
	fmt.Println("Sim::Simulacion")

	var pp, dttt float64
	ti := time.Now()

	// Ensure celdas is initialized once before the simulation loop
	//s.IniciarCel()

	for ic := 0; ic < s.nc; ic++ {
		s.icc = ic
		if ic%s.ncvc == 0 && s.metodo == "vc" {
			s.calcularVC()
		}

		for ia := 0; ia < len(s.atoms); ia++ {
			for id := 0; id < s.nd; id++ {
				pp = s.atoms[ia].P[id] + s.dt*s.atoms[ia].V[id] + s.dt*s.dt*s.atoms[ia].A[id]*0.5

				if pp > s.L {
					pp -= s.L
				}
				if pp < 0 {
					pp += s.L
				}
				s.atoms[ia].P[id] = pp
			}
		}

		for ia := 0; ia < len(s.atoms); ia++ {
			for id := 0; id < s.nd; id++ {
				s.atoms[ia].V[id] += s.dt * s.atoms[ia].A[id] / 2.0
			}
		}

		switch s.metodo {
		case "vc":
			s.AceleracionesVC()
		case "orig":
			s.Aceleraciones()
		case "cel":
			s.CalcularCel()
			s.AceleracionesCel()
		}

		for ia := 0; ia < len(s.atoms); ia++ {
			for id := 0; id < s.nd; id++ {
				s.atoms[ia].V[id] += s.dt * s.atoms[ia].A[id] / 2.0
			}
		}

		if ic == 0 {
			fmt.Println("ic, temp, dens, et, ec, ei, dttt")
		}

		if ic%s.nci == 0 {
			s.Prop()
			tf := time.Now()
			dttt = tf.Sub(ti).Seconds()
			fmt.Printf("%d %f %f %f %f %f %f\n", ic, s.temp, s.dens, s.et, s.ec, s.ei, dttt)
		}
	}
}

func main() {
	rand.Seed(time.Now().UnixNano())
	sim := Sim{}
	sim.Datos()
	fmt.Println("Los datos actuales son, Densidad:  ", sim.dens)
	sim.EscribirDatos()
	sim.IniciarAtomos()
	if sim.cuad == "cuad" {
		sim.Cuadrada()
	}
	// Inic
	//sim.Aceleraciones()
	switch sim.metodo {
	case "vc":
		sim.calcularVC()
		sim.AceleracionesVC()
	case "orig":
		sim.Aceleraciones()
	case "cel":
		sim.IniciarCel()
		sim.CalcularCel()
		sim.AceleracionesCel()
	}
	sim.Simulacion()
}
