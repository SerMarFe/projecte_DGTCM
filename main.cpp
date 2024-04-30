#include <iostream>
#include <cmath>
using namespace std;

//----------------------
// Entrades
//----------------------

// Físiques
const float e = 1; // [m] Gruix de la paret
const float S = 1; // [m^2] Secció de la paret
const float rho = 2400; // [kg/m^3] Densitat del fluid
const float lambda = 220; // [W/(m*K)] Conductivitat tèrmica de la paret (coeficient de transferència de calor per conducció)
const float alpha_ext = 8.7; // [W/(m^2*K)] Coeficient de transferència de calor per convecció exterior
const float T_0 = 21.0; // [ºC] Temperatura inicial de la paret
const float A_0 = 21.0; // [ºC] Temperatura inicial exterior
const float A_1 = 4.3; // [ºC] Amplitud de la variació de la temperatura segons l'hora del dia
const float A_2 = 7.5; // [ºC] Amplitud de la variació de la temperatura segons el dia de l'any
const float T_w = 21.0; // [ºC] Temperatura de la paret de la dreta
float T_ext; // [ºC] Temperatura exterior, funció del temps
const float w1 = 2*M_PI/(24*3600); // freqüència 1 d'oscil·lació de la temperatura
const float w2 = 2*M_PI/(365*24*3600); // freqüència 2 d'oscil·lació de la temperatura

// Numèriques
const int N = 10; // Nombre de volums de control de la discretització
const int dt = 3600; // [s] increment de temps de discretització temporal
const int t_fi = 100*24*3600; // [s] temps final
const int nT = 2; // Factor per reduir el nombre d'instants de temps en què guardem el mapa de temperatures
const float beta = 0.5; // Model --> 0: explícit, 0.5: Crank-Nicolson, 1: implícit
const float delta = 1e-3; // Criteri de convergència


//----------------------
// Definició de vectors
//----------------------

// Si tenim N volums de control tindrem N+1 nodes donada la geometria del problema
const float xcv[N] = {0}; // [m] Posició de les cares dels volums de control (l'array comença amb 0 per tant N --> N+1 elements)
const float xe[N] = {0}; // [m] Posició de les cares east
const float xw[N] = {0}; // [m] Posició de les cares west
const float xp[N+1] = {0}; // [m] Posició dels nodes
const float V = S*e/N; // [m^3] Volum dels volums de control
float ap[N+1] = {0}; // Coeficients de discretització dels nodes
float ae[N+1] = {0}; // Coeficients de discretització cara east
float aw[N+1] = {0}; // Coeficients de discretització cara west
float bp[N+1] = {0}; // Coeficients de discretització
float T[N+1][2] = {0}; // [ºC] Mapa de temperatures en instants contigus (n, n+1)
float T_hist[N+1][int((t_fi/dt)/nT)] = {0}; // [ºC] Mapa de temperatures a guardar


//----------------------
// Mapa inicial
//----------------------

void mapa_inicial() {
    for (int i = 0; i < N; i++) {
        T[i][1] = (T_0+T_w)/2;
    }
}


int main() {
    mapa_inicial();

    for (int i = 0; i<N; i++){
        cout<<T[i][1]<<endl;
    }
    return 0;
}