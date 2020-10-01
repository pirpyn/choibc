#include "../hoibc/hoibc.hpp"
#include "../hoibc/hoibc_bessel.hpp"
#include "../main/dump_csv.hpp"
#include "../main/dump_csv.cpp"
#include <iostream>

using namespace std;
using namespace hoibc;

int main(int argc, char* argv[]) {
  cout << "Testing value of impedance cylindric one layer in \"impedance_cylinder_verif.csv\"" << endl;

  hoibc::data_t data;
  data.main.frequency = 0.2;
  data.main.s1 = { 0., 1., 0.1};
  data.main.s2 = { 0., 0.9, 0.1};

  data.material.thickness = hoibc::array<hoibc::real>({ 0.05 });
  data.material.epsr      = hoibc::array<hoibc::complex>({ hoibc::complex(1.,-1.) });
  data.material.mur       = hoibc::array<hoibc::complex>({ hoibc::complex(1.,0.) });

  const hoibc::real k0 = hoibc::free_space_wavenumber(data.main.frequency);

  hoibc::hoibc_ibc3 ibc;
  ibc.mode = hoibc::mode_t::Z;
  ibc.suc = false;
  ibc.normalised = true;
  ibc.inner_radius = 1.;
  ibc.outer_radius = ibc.inner_radius + data.material.thickness[0];
  const hoibc::real r0 = ibc.inner_radius;
  const hoibc::real r1 = ibc.outer_radius;

  // Set up the scan of incident angle
  // For the plane, depends on (kx,ky)
  // For the cylinder, depends on kz, incidence is from theta = 0, but expressed as a Fourier serie over n, truncated
  // For the sphere, incidence is always from theta, phi = 0, but expressed as a Mie serie over n, truncated
  hoibc::array<hoibc::real> vn, vkz, s1 ,s2;
  ibc.set_fourier_variables(data,vn,vkz,s1,s2);

  hoibc::big_matrix<hoibc::complex> impedance_ex = hoibc::cylinder::impedance_infinite(vn,vkz,k0,data.material,r0);
  hoibc::big_matrix<hoibc::complex> impedance_ap = hoibc::big_init<hoibc::complex>(vn.size(),vkz.size(),hoibc::complex(0.,0.));

  const hoibc::complex eps = data.material.epsr[0];
  const hoibc::complex mu = data.material.mur[0];
  const hoibc::complex eta = std::sqrt(mu/eps);
  const hoibc::complex k = k0*std::sqrt(eps*mu);
  for (unsigned int i = 0; i < vn.size(); i++)
  {
    const hoibc::real n = vn[i];
    for  (unsigned int j = 0; j < vkz.size(); j++)
    {
      const hoibc::real kz = vkz[j];
      const hoibc::complex k3 = std::sqrt(k*k - kz*kz);

      hoibc::complex Sn =
        (hoibc::bessel2(n,k3*r1)*hoibc::bessel1p(n,k3*r0) - hoibc::bessel1(n,k3*r1)*hoibc::bessel2p(n,k3*r0))
        /
        (hoibc::bessel2p(n,k3*r1)*hoibc::bessel1p(n,k3*r0) - hoibc::bessel1p(n,k3*r1)*hoibc::bessel2p(n,k3*r0));

      hoibc::complex Tn =
        (hoibc::bessel2p(n,k3*r1)*hoibc::bessel1(n,k3*r0) - hoibc::bessel1p(n,k3*r1)*hoibc::bessel2(n,k3*r0))
        /
        (hoibc::bessel2(n,k3*r1)*hoibc::bessel1(n,k3*r0) - hoibc::bessel1(n,k3*r1)*hoibc::bessel2(n,k3*r0));

      impedance_ap[i][j][0][0] = ci*k3/(k*eta)*Sn;
      impedance_ap[i][j][1][0] = ci*kz*n/(k*k3*eta*r1)*Sn;
      impedance_ap[i][j][0][1] = ci*kz*n/(k*k3*eta*r1)*Sn;
      impedance_ap[i][j][1][1] = -ci*k/(k3*eta)*Tn + ci*n*n*kz*kz / (k*k3*k3*k3*eta*r1*r1)* Sn;
      impedance_ex[i][j] = inv(impedance_ex[i][j]);
    }
  }
  const std::string filename = "impedance_cylinder_verif.csv";
  // using hoibc::operator-;
  dump_to_csv(filename, vn, s2, impedance_ap - impedance_ex, "n", "kz/k0", "Z_O-Z_P", "");

  const bool check = abs((impedance_ap - impedance_ex)[0][0][0][0])<1e-8;

  if (!check) {
    cout << "[FAIL] " << "|(Z_O-Z_P)[0][0][0][0]| > 1e-8" << endl;
    ifstream myfile;
    myfile.open(filename);
    string fileline;
    for (int i=0; i < 2; i++) {
      getline(myfile,fileline);
      cout << fileline << endl;
    }

    myfile.close();
  } else {
    cout << "[OK?] " << "|(Z_O-Z_P)[0][0][0][0]| < 1e-8" << endl;
  }
  return (!check);
}