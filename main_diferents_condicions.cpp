#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

//----------------------
// Entrades
//----------------------

// Físiques
const double e = 1.5; // [m] Gruix de la paret
const double S = 1; // [m^2] Secció de la paret
const double rho = 760; // [kg/m^3] Densitat del fusta de roure
const double cp = 2384.88; // [J/(kg*K)] Capacitat calorífica a pressió constant del fusta de roure
const double lambda = 0.16; // [W/(m*K)] Conductivitat tèrmica de la fusta de roure (coeficient de transferència de calor per conducció)
const double alpha_ext = 8.7; // [W/(m^2*K)] Coeficient de transferència de calor per convecció exterior
const double T_0 = 21.0; // [ºC] Temperatura inicial de la paret
const double A_0 = -15.0; // [ºC] Temperatura inicial exterior
const double A_1 = 12.0; // [ºC] Amplitud de la variació de la temperatura segons l'hora del dia
const double A_2 = 20.0; // [ºC] Amplitud de la variació de la temperatura segons el dia de l'any
const double T_w = 21.0; // [ºC] Temperatura de la paret de la dreta
double T_ext[2] = {}; // [ºC] Temperatura exterior, funció del temps, definim dues posicions, per al temps actual i per a l'anterior
const double w1 = 2*M_PI/(24*3600); // freqüència 1 d'oscil·lació de la temperatura
const double w2 = 2*M_PI/(365*24*3600); // freqüència 2 d'oscil·lació de la temperatura

// Numèriques
const int N = 51; // Nombre de volums de control de la discretització
const int Nnodes = N+2; // Nombre de nodes
const int Ncares = N+1; // Nombre de cares 
const double dt = 10; // [s] increment de temps de discretització temporal
const int t_fi = 100*24*3600; // [s] temps final
const int freq = 1; // Factor per reduir el nombre d'instants de temps en què guardem el mapa de temperatures. Corresponent a: cada quants instantst guardem? Si freq=1 guardem tot; si 2, la meitat...
const double Beta = 0.5; // Model --> 0: explícit, 0.5: Crank-Nicolson, 1: implícit
const double fr = 1; // Factor de relaxació;
const double delta = 1e-12; // Criteri de convergència

//----------------------
// Definició de vectors
//----------------------

// Si tenim N volums de control tindrem N+2 nodes donada la geometria del problema
double dx = e/N; // [m] Distància entre nodes
double xcv[Ncares] = {}; // [m] Posició de les cares dels volums de control
double xe[Ncares] = {}; // [m] Posició de les cares east
double xw[Ncares] = {}; // [m] Posició de les cares west
double xp[Nnodes] = {}; // [m] Posició dels nodes
double dpw[Nnodes] = {}; // [m] Distància entre cares WP
double dpe[Nnodes] = {}; // [m] Distància entre cares PE
double V = S*dx; // [m^3] Volum dels volums de control
double ap[Nnodes] = {}; // Coeficients de discretització dels nodes
double ae[Nnodes] = {}; // Coeficients de discretització cara east
double aw[Nnodes] = {}; // Coeficients de discretització cara west
double bp[Nnodes] = {}; // Coeficients de discretització
double T[Nnodes][2] = {}; // [ºC] Mapa de temperatures en instants contigus (n, n+1)
double T_est[Nnodes][2] = {}; // [ºC] mapa de temperatures estimat
double T_hist[Nnodes][1+int((t_fi/dt)/freq)] = {}; // [ºC] Mapa de temperatures a guardar

// Variables de comprovació i calculs
double Qin[int(t_fi/dt)/freq] = {};
double Qout[int(t_fi/dt)/freq] = {};
double Qacc[int(t_fi/dt)/freq] = {};
double Qtotal[int(t_fi/dt)/freq] = {};
double T_analitica[Nnodes] = {}; // Valor analític de la temperatura en règim permanent si temperatura externa constant
//------------------------------------------
// Definicó dels prototipus de les funcions
//------------------------------------------

void posicions_nodes();
void mapa_inicial();
void calcula_coefs_disc();
bool comprova_convergencia();
void actualitza_T_est(int j, int k);
void actualitza_T();
void gauss_seidel();
void nou_temps();
void calcula_T_ext(int i);

void balanc_global(int i);
void solucio_analitica();
void exporta_resultats();

//----------------------
// Main
//----------------------

int main() {
    posicions_nodes();
    mapa_inicial();
    nou_temps();
    cout<< setprecision(12);
    cout<<"Temperatura dels nodes:"<<endl;
    for (int i = 0; i<Nnodes; i++){
        cout<<"Node "<<i<<" (x = "<<xp[i]<<"m): "<<T[i][1]<<"ºC"<<endl;
    };
    solucio_analitica();
    exporta_resultats();
    return 0;
}

//-----------------------
// Definició de funcions
//-----------------------

void posicions_nodes() {
    xcv[0] = 0;
    xe[0] = xcv[0];
    xw[0] = 0;
    xp[0] = 0;
    for (int i = 1; i < Ncares; i++) {
        xcv[i] = i*dx;
        xw[i] = xcv[i-1];
        xe[i] = xcv[i];
        xp[i] = (xw[i]+xe[i])/2;
    };
    xp[N+1] = e; 
    for (int i = 0; i < Ncares; i++) {
        dpe[i] = xp[i+1]-xp[i];
        dpw[i] = xp[i]-xp[i-1];
    };
    dpw[0] = 0; // no es fa servir (ja que w del primer node de la paret no existeix)
    dpw[N+1] = xp[N+1] - xp[N];
    // la posició dpe[N+1] no es fa servir, per això també val 0
};

void mapa_inicial() {
    for (int i = 0; i < Nnodes; i++) {
        T[i][0] = T_0; // Dada enunciat
        T_hist[i][0] = T_0;
    };
};

bool comprova_convergencia() {
    for (int i = 0; i < Nnodes; i++) {
        double dif = abs(T[i][1] - T_est[i][1]);
        if (dif > delta) return false;
    };
    return true;
};

void actualitza_T_est(int j, int k) {
    for (int i = 0; i < Nnodes; i++) {
        T_est[i][j] = T[i][k];
    };
};

void actualitza_T() {
    for (int i = 0; i < Nnodes; i++) {
        T[i][0] = T[i][1];
    };
};

void guarda_T(int j) {
    for (int i = 0; i < Nnodes; i++) {
        T_hist[i][j] = T[i][1];
    }
    
};

void gauss_seidel() {
    double conv = false;
    while (!conv) {
        //calcula_coefs_disc(); // Si variessin amb la temperatura caldria descomentar i adaptar-ho
        for (int i = 0; i < Nnodes; i++) {
            T[i][1] = (aw[i]*T_est[i-1][1] + ae[i]*T_est[i+1][1] + bp[i])/ap[i];
            T[i][1] = T[i][0] + fr*(T[i][1]-T[i][0]);
            //cout << "temperatura " << i << ": " << T[i][1] <<endl;
        };
        conv = comprova_convergencia();
        if (!conv) actualitza_T_est(1,1);
        //cout << "Condició: " << conv << endl;
    };
};

void nou_temps() {
    for (int i = 0; i < t_fi/dt; i++) {
        actualitza_T_est(1,0);
        calcula_T_ext(i);
        calcula_coefs_disc();
        gauss_seidel();
        if (i%freq==0) guarda_T(int(1+i/freq)); // Per guardar les dades cada 'freq' vegades increments de temps
        balanc_global(i);
        actualitza_T();
        cout<<i/(t_fi/dt)*100<<"%"<<endl;
    };
};

void calcula_coefs_disc() {
    aw[0] = 0;
    ae[0] = lambda*S/dpe[0];
    ap[0] = ae[0] + alpha_ext*S;
    bp[0] = alpha_ext*T_ext[1]*S;

    for (int i = 1; i < N+1; i++) {
        aw[i] = Beta*lambda*S/dpw[i];
        ae[i] = Beta*lambda*S/dpe[i];
        ap[i] = aw[i] + ae[i] + rho*V*cp/dt;
        bp[i] = rho*V*cp*T[i][0]/dt + (1-Beta)*(-lambda*(T[i][0]-T[i-1][0])*S/dpw[i] + lambda*(T[i+1][0]-T[i][0])*S/dpe[i]);
    };

    aw[N+1] = 0;
    ae[N+1] = 0;
    ap[N+1] = 1;
    bp[N+1] = T_w;
};

void calcula_T_ext(int i) {
    T_ext[0] = A_0 + A_1*sin(w1*i*dt) + A_2*sin(w2*i*dt);
    T_ext[1] = A_0 + A_1*sin(w1*(i+1)*dt) + A_2*sin(w2*(i+1)*dt); 
    //T_ext[0] = 700;
    //T_ext[1] = 700;
};

void balanc_global(int i) {
    // Qin = Potència de calor entrant per conducció a través de la paret de l'esquerra
    // Qout = Potència de calor de conducció de l'últim volum de control a l'últim node, corresponent a la superfície de la dreta
    // Qacc = Potència de calor acumulada pels volums de control
    // Com que es fa el balanç d'energia per unitat de temps, és a dir, per a cada step de temps, es fa el balanç de potència, cal usar l'esquema d'integració emprat mitjançant Beta.
    Qin[i] = -Beta*lambda*(T[1][0]-T[0][0])*S/dpe[0] + -(1-Beta)*lambda*(T[1][1]-T[0][1])*S/dpe[0];
    Qout[i] = -Beta*lambda*(T[N+1][0]-T[N][0])*S/dpe[N] + -(1-Beta)*lambda*(T[N+1][1]-T[N][1])*S/dpe[N];
    Qacc[i] = 0;
    for (int j=1; j < N+1; j++) {
        Qacc[i] = Qacc[i] + rho*V*cp*(T[j][1]-T[j][0])/dt;
    };
    Qtotal[i] = Qin[i] - Qout[i] - Qacc[i]; // Ha de valer 0
};

void solucio_analitica() {
    for (int i=0; i < Nnodes; i++)
    T_analitica[i] = (e*T_ext[1] + lambda*T_w/alpha_ext + (T_w-T_ext[1])*xp[i]) / (e+lambda/alpha_ext);
};

void exporta_resultats() {
    ofstream resultats("resultats/resultats_temperatura.csv");
    resultats << setprecision(8);
    if (!resultats.is_open()) {
        cerr << "No s'ha pogut obrir el fitxer per escriure-hi" << endl;
    };
    // Capçal
    resultats << "Temps final: " << t_fi << "\n";
    // Escriu la temperatura
    for (int i = 0; i < Nnodes; i++) {
        //resultats << i; // Index del node
        resultats << xp[i];
        for (int j = 0; j < int((t_fi/dt)/freq); j++) {
            resultats <<  "," << T_hist[i][j];
        };
        resultats << "\n";
    };

    ofstream resultats_analitics("resultats/resultats_analitics_temperatura.csv");
    if (!resultats_analitics.is_open()) {
        cerr << "No s'ha pogut obrir el fitxer per escriure-hi" << endl;
    };
    //Capçal
    resultats_analitics << "Temps final: " << t_fi << "\n";
    for (int i = 0; i < Nnodes; i++) {
        resultats_analitics << T_analitica[i] << ",";
    };

    ofstream balanc_calor("resultats/balanc_calor.csv");
    if (!balanc_calor.is_open()) {
        cerr << "No s'ha pogut obrir el fitxer per escriure-hi" << endl;
    };
    //Capçal
    balanc_calor << "Temps final: " << t_fi << "\n";
    balanc_calor << "t, Qin, Qout, Qacc, Qtotal" << "\n";
    for (int i = 0; i < int((t_fi/dt)/freq); i++) {
        balanc_calor << i*dt*freq << ",";
    };
    balanc_calor << "\n";    
    for (int i = 0; i < int((t_fi/dt)/freq); i++) {
        balanc_calor << Qin[i] << ",";
    };
    balanc_calor << "\n";
    for (int i = 0; i < int((t_fi/dt)/freq); i++) {
        balanc_calor << Qout[i] << ",";
    };
    balanc_calor << "\n";
    for (int i = 0; i < int((t_fi/dt)/freq); i++) {
        balanc_calor << Qacc[i] << ",";
    };
    balanc_calor << "\n";
    for (int i = 0; i < int((t_fi/dt)/freq); i++) {
        balanc_calor << Qtotal[i] << ",";
    };
};