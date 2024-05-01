#include <iostream>
#include <iomanip>
using namespace std;

//----------------------
// Entrades
//----------------------

// Físiques
const double e = 1; // [m] Gruix de la paret
const double S = 1; // [m^2] Secció de la paret
const double rho = 2400; // [kg/m^3] Densitat del fluid
const double lambda = 220; // [W/(m*K)] Conductivitat tèrmica de la paret (coeficient de transferència de calor per conducció)
const double alpha_ext = 8.7; // [W/(m^2*K)] Coeficient de transferència de calor per convecció exterior
const double T_w = 21.0; // [ºC] Temperatura de la paret de la dreta
double T_ext = 750; // [ºC] Temperatura exterior
// Numèriques
const int N = 50; // Nombre de volums de control de la discretització
const int Nnodes = N+2; // Nombre de nodes
const int Ncares = N+1; // Nombre de cares 
const double delta = 1e-14; // Criteri de convergència


//----------------------
// Definició de vectors
//----------------------

// Si tenim N volums de control tindrem N+1 cares i N+2 nodes donada la geometria del problema
double dx = e/N; // [m] Distància entre nodes
double xcv[Ncares] = {}; // [m] Posició de les cares dels volums de control
double xe[Ncares] = {}; // [m] Posició de les cares east
double xw[Ncares] = {}; // [m] Posició de les cares west
double xp[Nnodes] = {}; // [m] Posició dels nodes
double dpw[Nnodes] = {}; // [m] Distància entre cares WP
double dpe[Nnodes] = {}; // [m] Distància entre cares PE
double V = S*dx; // [m^3] Volum dels volums de control
double T[Nnodes]; // [ºC] mapa de temperatures
double T_est[Nnodes]; // [ºC] mapa de temperatures estimat
struct coefsDisc {
    double ap[Nnodes]; // Coeficients de discretització dels nodes
    double ae[Nnodes]; // Coeficients de discretització cara east
    double aw[Nnodes]; // Coeficients de discretització cara west
    double bp[Nnodes]; // Coeficients de discretització terme independent
};
coefsDisc a;

const double T0_analitica = (alpha_ext*T_ext+lambda*T_w/e)/(alpha_ext+lambda/e);

//------------------------------------------
// Definicó dels prototipus de les funcions
//------------------------------------------

void posicions_nodes();
void mapa_inicial();
void calcula_coefs_disc();
bool comprova_convergencia();
void actualitza_T_est();
void gauss_seidel();


//----------------------
// Main
//----------------------

int main() {
    posicions_nodes();
    mapa_inicial();
    calcula_coefs_disc();
    gauss_seidel();
    cout<< setprecision(12);
    cout<<"Resultat analític de la temperatura del primer node: "<<endl;
    cout<<"T0 = "<<T0_analitica<<" ºC (analític)"<<endl;
    cout<<endl;
    cout<<"Temperatura dels nodes:"<<endl;
    for (int i = 0; i<Nnodes; i++){
        cout<<"Node "<<i<<" (x = "<<xp[i]<<"m): "<<T[i]<<"ºC"<<endl;
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
    dpw[N+1] = xp[N+1] - xp[N];
};

void mapa_inicial() {
    for (int i = 0; i < Nnodes; i++) {
        T_est[i] = (T_ext+T_w)/2; // Com a primera aproximació, inicialitzem el mapa de temperatures al valor promig
    };
};

bool comprova_convergencia() {
    for (int i = 0; i < Nnodes; i++) {
        double dif = abs(T[i]-T_est[i]);
        if (dif > delta) return false;
    };
    return true;
};

void actualitza_T_est() {
    for (int i = 0; i < Nnodes; i++) {
        T_est[i] = T[i];
    };
};

void gauss_seidel() {
    double conv = false;
    while (!conv) {
        for (int i = 0; i < Nnodes; i++) {
            T[i] = (a.aw[i]*T[i-1] + a.ae[i]*T_est[i+1] + a.bp[i])/a.ap[i]; //a.aw[i]*T[i-1] o a.aw[i]*T_est[i-1] ?? (funciona igual però)
            //cout << "temperatura " << i << ": " << T[i] <<endl;
        };
        conv = comprova_convergencia();
        if (!conv) actualitza_T_est();
        //cout << "Condició: " << conv << endl;
    };
};

void calcula_coefs_disc() {
    a.aw[0] = 0;
    a.ae[0] = lambda*S/dpe[0];
    a.ap[0] = a.ae[0] + alpha_ext*S;
    a.bp[0] = alpha_ext*T_ext*S;

    for (int i = 1; i < N+1; i++) {
        a.aw[i] = lambda*S/dpw[i];
        a.ae[i] = lambda*S/dpe[i];
        a.ap[i] = a.aw[i] + a.ae[i];
        a.bp[i] = 0;
    };

    a.aw[N+1] = 0;
    a.ae[N+1] = 0;
    a.ap[N+1] = 1;
    a.bp[N+1] = T_w;
};