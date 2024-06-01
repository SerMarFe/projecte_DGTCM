#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

//----------------------
// Entrades
//----------------------

// Físiques
const double e = 1; // [m] Gruix de la paret
const double S = 1; // [m^2] Secció de la paret
const double rho = 2400; // [kg/m^3] Densitat del fluid
const double cp = 900; // [J/(kg*K)] Capacitat calorífica a pressió constant
const double lambda = 220; // [W/(m*K)] Conductivitat tèrmica de la paret (coeficient de transferència de calor per conducció)
const double alpha_ext = 8.7; // [W/(m^2*K)] Coeficient de transferència de calor per convecció exterior
const double T_0 = 21.0; // [ºC] Temperatura inicial de la paret
const double A_0 = 21.0; // [ºC] Temperatura inicial exterior
const double A_1 = 4.3; // [ºC] Amplitud de la variació de la temperatura segons l'hora del dia
const double A_2 = 7.5; // [ºC] Amplitud de la variació de la temperatura segons el dia de l'any
const double T_w = 21.0; // [ºC] Temperatura de la paret de la dreta
double T_ext[2] = {}; // [ºC] Temperatura exterior, funció del temps, definim dues posicions, per al temps actual i per a l'anterior
const double w1 = 2*M_PI/(24*3600); // freqüència 1 d'oscil·lació de la temperatura
const double w2 = 2*M_PI/(365*24*3600); // freqüència 2 d'oscil·lació de la temperatura

// Numèriques
const int N = 5; // Nombre de volums de control de la discretització
const int Nnodes = N+2; // Nombre de nodes
const int Ncares = N+1; // Nombre de cares 
const double dt = 1; // [s] increment de temps de discretització temporal
const int t_fi = 1000;//100*24*3600l; // [s] temps final
const int freq = 1; // Factor per reduir el nombre d'instants de temps en què guardem el mapa de temperatures. Corresponent a: cada quants instantst guardem? Si freq=1 guardem tot; si 2, la meitat...
const double Beta = 0.5; // Model --> 0: explícit, 0.5: Crank-Nicolson, 1: implícit
const double fr = 100; // Factor de relaxació;
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
double T_est[Nnodes][2] = {}; // [ºC] mapa de temperatures estimat //només fem servir la posició 2 corresponent a l'instant n+1 (és a dir, T_est[i][0] no es fa servir mai)
double T_hist[Nnodes][int((t_fi/dt)/freq)] = {}; // [ºC] Mapa de temperatures a guardar

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
        if (i%freq==0) guarda_T(1+i/freq); // Per guardar les dades cada 'freq' vegades increments de temps
        actualitza_T();
        cout<<i/(t_fi/dt)*100<<"%"<<endl;
    };
};

void calcula_coefs_disc() {
    aw[0] = 0;
    ae[0] = Beta*lambda*S/dpe[0];
    ap[0] = ae[0] + Beta*alpha_ext*S;
    bp[0] = Beta*alpha_ext*T_ext[1]*S + (1-Beta)*(alpha_ext*(T_ext[0]-T[0][0])*S + lambda*(T[1][0]-T[0][0])*S/dpe[0]);

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
    //T_ext[0] = A_0 + A_1*sin(w1*i*dt) + A_2*sin(w2*i*dt);
    //T_ext[1] = A_0 + A_1*sin(w1*(i+1)*dt) + A_2*sin(w2*(i+1)*dt); 
    T_ext[0]=750;
    T_ext[1]=750;
};

// REVISAR COEFS DISC. Algo falla, mirar el valor de T node 0 a cada iteració, és estrany.