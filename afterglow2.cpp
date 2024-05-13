// Primary changes to this code are the smooth function inclusions. This code does not contain the x-ray integral
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::setw;
using std::cerr;
using std::max;

double Ej =5.0E+51; // Initial energy in ergs
double n = 0.05; // particle density
double e = 4.8032E-10; //charge of electron
double ep_e = 0.01; // E- field dispersion
double ep_B = pow(10,-3.0); // B-field Dispersion
double p = 2.17; // power index
double t = 30.0; // time
double G1 = 100; // Initial gamma from prompt emission
double m_p = 1.67E-27; // mass of proton
double m_e = 9.109E-28; // mass of electron
double c = 29979245800; // speed of light in cm/s
double t_sec = 6.65E-25; //thompson cross section
double t1 = .1;
double t2 = .1;
double delta_x = 1.0;
double delta_y = 1.0;
double sum = 0.0;
double D = 1.2E+26; // distance from black hole
double nu = 3.0E+9; // radio frequency GHz
double xr = 2.417990504024e+17; //1KeV Ghz
double g = 0.0; // initialize gamma
double b = 0.0; // intialize beta
double sum_fx = 0.0;
double total_fx = 0.0;
double h_max = 20.0; // max theta
double p_max = 12.0; // max phi

double rdeg(double x) // Convert to radians & initialize static T_obs
{
    return (M_PI / 180.0) * x;
}
double x = rdeg(20.0); // .34906
double y = rdeg(25.0); // .52359
double z = rdeg(12); //  .20944

double lcos(double x,double y,double z) // observer angle to be set to y = 30 off-axis observer equation inputs are constant for now

{
	return cos(x) * cos(y) + sin(x) * sin(y) * cos(z);
}


double G(double x) // Initial Gamma: takes in an angle in radians
{
  return G1*(0.148*pow(x,-1*0.15364) - 0.01781);
  //return G1*(0.67148818*(exp(1/(12.258359390)*-x)) + 1.02383338);
}

double B(double g) //  Initial Beta - Gamma dependent relation takes in G(x):g
{
    return pow((1.0 - (1.0 / pow(g, 2.0))), 0.5);    
}


double E(double x) //Energy curve takes in radians
{
  return Ej*(0.766989*pow(x,-1*0.15341) + 0.93639);
  //return Ej*(0.12114057*(exp(1/(9.745295)*-x)) + 0.000541651);
}

double Bp(double g) // B-prime takes in functional input is the function G(t):g
{
	return pow(32*M_PI*ep_B*g*(g - 1.0)* n *m_p*pow(c,2.0),0.5);	
}

double Ga( double b) // Gamma actual - strictly for Gamma solution takes in B(x):b
{
	return 1.0/(pow(1.0 - pow(b,2.0),0.5));
}

double R_dcel(double x) // function  theta E0 is piece-wise defined for regions greater than t1 = .1
{
	return pow((3.0 * E(x)) / (4.0 * M_PI * n * m_p * pow(c, 2.0) * pow(G(x), 2.0) * pow(B(G(x)), 2.0)), 0.333);
}

double B_dcel(double x, double r) // Deceleration Takes in an angle for x - this is constant. Takes in a radius: r calculated numerically.
{	
	if (R_dcel(x) >= r){ 
	return B(G(x));
	}
	else {
	return (pow((R_dcel(x)/r),1.5)*pow(G(x),1.0)*pow(B(G(x)),1.0))/(pow(1.0 + pow(R_dcel(x)/r,3.0)*pow(G(x),2.0)*pow(B(G(x)),2.0),0.5));
	}
}

double g_m(double g) // Gamma-max takes G(t):g as an input
{
	if( ((p-2.0)/(p-1.0))*ep_e*(g - 1.0)*(m_p/m_e) >= 1){
		return ((p-2.0)/(p-1.0))*ep_e*(g - 1.0)*(m_p/m_e);
	}
	else{
		{
	return 1;
	}
}
}
double del(double g, double b, double x, double y, double z) // Doppler Factor - Takes Gamma-0 and Beta-0: g and b as a functional input
{
	return 1.0/(g*(1.0 - b*lcos(x,y,z)));
}

double nu_p_m (double b, double g) // Synchtron Max takes in Bp:b and g_m:g as functional inputs
{
	return (3.0*e*b*pow(g,2.0))/(4.0*M_PI*m_e*c);	
}

double nu_p(double d) // Spectral power up to nu'syn takes in gamma and beta as well as the constant angle. Takes in functional del:d chekc gamma-m < unity
{
	return nu/d;
}
double xr_p(double d) // Spectral power up to nu'syn takes in gamma and beta as well as the constant angle. Takes in functional del:d chekc gamma-m < unity
{
	return xr/d; // x-ray input
}

double p_p_nu_p_m(double b) // P'nu'max takes in Bp:b as functional inputs
{
	return (2.0*m_e*t_sec*pow(c,2.0)/(9.0*e))*b;
}

double p_p_nu_p(double g, double b,double d) // P'nu' takes in Gamma, Beta, and Delta as functional inputs: g, b, and d respectively. 
{
    double nu_p_val = nu_p(d);  // Calculate nu_p value
    
    if (nu_p_val >= nu_p_m(b,g)) {
        return p_p_nu_p_m(b) * pow(nu_p_val / nu_p_m(b,g), (1.0 - p) / 2.0);
    } else {
        return p_p_nu_p_m(b) * pow(nu_p_val / nu_p_m(b,g), 0.333);
    }
}
double p_p_x_p(double g, double b,double d) // P'nu' takes in Gamma, Beta, and Delta as functional inputs: g, b, and d respectively. 
{
    double x_p_val = xr_p(d);  // Calculate nu_p value
    
    if (x_p_val >= nu_p_m(b,g)) {
        return p_p_nu_p_m(b) * pow(x_p_val / nu_p_m(b,g), (1.0 - p) / 2.0);
    } else {
        return p_p_nu_p_m(b) * pow(x_p_val / nu_p_m(b,g), 0.333);
    }
}

double func(double t, double r) // dr/dt = Bc/1-Blcos
{
    double bdcel = B_dcel(t, r);
    return (bdcel*c) / (1.0 - bdcel * lcos(x, y, z));
}

int main()
{
	
    double minT = 0.01; // test values
    double maxT = 1000.0;
    double dt = 1.0;
    double r0 = pow(10,15); 
    int steps = std::lround((maxT-minT)/dt+1);
    
    
    std::vector<double> tvals(steps);
	std::vector<double> bofg(steps);
    std::vector<double> roft(steps);
	std::vector<double> sums(steps);
	std::vector<double> bofr(steps);
	std::vector<double> goft(steps);
	std::vector<double> doft(steps);
	std::vector<double> Bpoft(steps);
	std::vector<double> nupoft(steps);
	std::vector<double> gmoft(steps);
	std::vector<double> nupmoft(steps);
	std::vector<double> xrpmoft(steps);
	std::vector<double> ppnumpoft(steps);
	std::vector<double> pxnumpoft(steps);
	std::vector<double> ppnup(steps);
	std::vector<double> ppxp(steps);
	std::vector<double> fx(steps);
	std::vector<double> fxr(steps);
	std::vector<double>sum_fx(steps,0.0);
	std::vector<double>sum_fxr(steps, 0.0);
	
	
    std::ofstream sol("roft2.txt");
    if(sol.is_open())
        cout << "sol stream opened successfully." << endl;
    else
    {
        cout << "sol stream failed to open!" << endl;
        return EXIT_FAILURE;
    }
	
    for (double x = rdeg(0.0) ; x <= h_max; x += delta_x) { //theta
			for (double z = 0.0; z <= 1; z += delta_y) {//phi
			
    tvals.at(0) = minT;
    roft.at(0) = r0; //initiallized for roft
	bofr.at(0) = B_dcel(rdeg(x),r0); // take from B_dcel
	goft.at(0) = Ga(bofr.at(0));
	doft.at(0) = del(goft.at(0),bofr.at(0),rdeg(x),y,rdeg(z));
	Bpoft.at(0) = Bp(goft.at(0));
	nupmoft.at(0) = nu_p_m(Bpoft.at(0),gmoft.at(0));
	xrpmoft.at(0) = nu_p_m(Bpoft.at(0),gmoft.at(0));
	ppnumpoft.at(0) = p_p_nu_p_m(Bpoft.at(0));
	fx.at(0) = 1.0*pow(10,23)*(1 / (M_PI*pow(D,2.0))*goft.at(0)*(pow(roft.at(0), 3.0) / 3.0)) * pow(doft.at(0),3)*ppnup.at(0)*sin(rdeg(z))*rdeg(x)*rdeg(z);
	fxr.at(0) = 1.0*pow(10,23)*(1 / (M_PI*pow(D,2.0))*goft.at(0)*(pow(roft.at(0), 3.0) / 3.0)) * pow(doft.at(0),3)*ppxp.at(0)*sin(rdeg(z))*rdeg(x)*rdeg(z);
    
    sol << tvals.at(0) << "\t" << roft.at(0) << "\t" << bofr.at(0) << "\t" << goft.at(0) << "\t" << doft.at(0) << "\t" << Bpoft.at(0) << "\t" << nupoft.at(0) << "\t" << ppnumpoft.at(0) << "\t" << fx.at(0) << endl;
   // int end_patch = fx.size();
	
    for(int i=1; i<steps; i++)//change deltat for log
    {
		
        tvals.at(i) = tvals.at(i-1)+dt;
        roft.at(i) = roft.at(i-1) + (3600.0*24.0)*dt*func(rdeg(x),roft.at(i-1));
		
		
		//for (int j = 0; j <= end_patch - 1; j++){
		//	fx.at(j) = 0;
		//}
		
		//for (double x = 0.0; x <= 20; x += delta_x) { 
		//	for (double y = 0.0; y <= 30; y += delta_y) {
				
		bofr.at(i) = B_dcel(rdeg(x), roft.at(i));
		goft.at(i) = Ga(bofr.at(i));
		doft.at(i) = del(goft.at(i),bofr.at(i),rdeg(x),y,rdeg(z)); // Check shift later
		Bpoft.at(i) = Bp(goft.at(i));	
		nupoft.at(i) = nu_p(doft.at(i)); 
		gmoft.at(i) = g_m(goft.at(i));
		nupmoft.at(i) = nu_p_m(Bpoft.at(i),gmoft.at(i));
		xrpmoft.at(i) = nu_p_m(Bpoft.at(i),gmoft.at(i));
		ppnumpoft.at(i) = p_p_nu_p_m(Bpoft.at(i));
		ppnup.at(i) = p_p_nu_p(Bpoft.at(i),goft.at(i),doft.at(i));
		ppxp.at(i) = p_p_x_p(Bpoft.at(i),goft.at(i),doft.at(i));
		  // }
		//}
		fxr.at(i) = (1/2.0E-20)*(1 / (4*M_PI*pow(D,2.0))*n*(pow(roft.at(i), 3.0) / 3.0)) * pow(doft.at(i),3)*ppxp.at(i)*sin(rdeg(z))*rdeg(z)*rdeg(x);//theta-phi
        fx.at(i) = (1/2.0E-27)*(1 / (4*M_PI*pow(D,2.0))*n*(pow(roft.at(i), 3.0) / 3.0)) * pow(doft.at(i),3)*ppnup.at(i)*sin(rdeg(z))*rdeg(z)*rdeg(x);//theta-phi
		sum_fx[i] += fx.at(i);
		sum_fxr[i] += fxr.at(i);
		
		
		//sol << tvals.at(i);
		//for (int j = 0;  j <= end_patch-1; j++)
		//	{ sol << "\t" << fx.at(j);
		//	}
		//sol << endl;
        //dt *= 1;
      

        sums.at(i) = sum;
       sol << tvals.at(i) << "\t" << roft.at(i) <<  "\t" << gmoft.at(i) << "\t" << goft.at(i) << "\t" << doft.at(i) << "\t" << Bpoft.at(i)  << "\t" << sum_fxr[i] << "\t" << sum_fx[i] << "\t" << fx.at(i) << endl;
	sum = 0.0;
	
	//cout << "rdeg" << lcos(x,y,z) << endl;
    //cout << "tvals[" << i << "]: " << tvals.at(i) << endl;
    //cout << "roft[" << i << "]: " << roft.at(i) << endl;
	//cout << "R_dcel[" << i << "]: " << R_dcel(x) << endl;
	//cout << "B_dcel[" << i <<"]: " << bofr.at(i) << endl;
    //cout << "goft[" << i << "]:" << goft.at(i)<< endl;
	//cout << "Bpoft[" << i << "]:" << Bpoft.at(i)<< endl;
	//cout << "nupmoft[" << i << "]:" << nupoft.at(i)<< endl;
	//cout << "gmoft[" << i << "]:" << gmoft.at(i)<< endl;
	//cout << "nupmoft[" << i << "]:" << nupmoft.at(i)<< endl;
	//cout << "ppnumpoft[" << i << "]:" << ppnumpoft.at(i)<< endl;
	//cout << "ppnup[" << i << "]:" << ppnup.at(i)<< endl;
	//cout << "fx[" << i << "]: " << fx.at(i) << endl;
	//cout << "Sum Final" << ":" << sum_fx[i] << endl;
		   }
		}
	}
	
	cout << "Data has been generated successfully." << endl;
    sol.close();
	
	ofstream last_entries("sum_fx2.txt");
    if (last_entries.is_open()) {
        // Starting index
        int start_index = max(0, static_cast<int>(sum_fx.size()) - 1000);
        // Write last entries only
        for (double i = start_index; i < sum_fx.size(); ++i) {
            last_entries<< tvals.at(i) << "\t" << sum_fx[i] << "\t" << sum_fxr[i] << "\n";
        }
        
       
        last_entries.close();
        cout << "Sum entries written" << endl;
    } else {
        cerr << "Unable to open file for writing!" << endl;
    }


    return EXIT_SUCCESS;
}

