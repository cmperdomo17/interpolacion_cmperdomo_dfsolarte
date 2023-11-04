#ifndef SPLINE3_H
#define SPLINE3_H

#include "util.h"

#include <cmath>
#include <vector>

using std::vector;
using util::gauss;

namespace interpolacion {
    class spline3 {
        public: 
            spline3(vector <double> p_x,
                    vector <double> p_y): x(p_x), y(p_y){
                    // Calcular las segundas derivadas
                    f2 = calcular_f2(); 
            }
            
            /**
             * @brief Evaluar el polinomio del trazador cúbico en x_int
             * @param x_int Punto a evaluar
             * @return Valor interpolado en x_int
            */
            double interpolar(double x_int){
                
                int i = 0; 
                int n = x.size(); /*!< Numero de datos*/
                int intervalos = n - 1; /*!< Numero de intervalos*/

                // Verificar que x_int esté dentro del rango de los datos
                if (x_int < x[0]) {
                    // x_int está por debajo del rango de los datos
                    return NAN;
                }
                if (x_int > x[n-1]) {
                    // x_int está por encima del rango de los datos
                    return NAN;
                }

                // Determinar el intervalo i en donde se encuentra x_int
                for (i = 1; i < intervalos; i++){
                    if (x_int >= x[i - 1] && x_int <= x[i]){
                        break;
                    }
                }

                // Evaluar el polinomio del trazador en x_int (18.36)

                double h = x[i] - x[i - 1];

                double a = ((f2[i - 1] / (6.0f * h)) * (pow(x[i] - x_int, 3.0f))) + ((f2[i] / (6.0f * h)) * (pow(x_int - x[i - 1], 3.0f)));

                double b = ((y[i - 1] / h) - ((f2[i - 1] * h) / 6.0f)) * (x[i] - x_int);

                double c = ((y[i] / h) - ((f2[i] * h) / 6.0f)) * (x_int - x[i - 1]);

                double resultado = a + b + c;

                return resultado;

            }

            /**
             * @brief Interpolar los coeficientes del trazador cúbico en x_int y mostrar el polinomio de cada subintervalo
             * @param x_int Punto a evaluar
             * @return Vector de coeficientes del polinomio
            */
            vector <double> interpolar_trazador(double x_int){
                
                int i = 0;
                int n = x.size(); /*!< Numero de datos*/
                int intervalos = n - 1; /*!< Numero de intervalos*/

                // Determinar el intervalo i en donde se encuentra x_int
                for (i = 1; i < intervalos; i++){
                    if (x_int >= x[i - 1] && x_int < x[i]){
                        break;
                    }
                }

                // Evaluar el polinomio del trazador en x_int (18.36)

                double h = x[i] - x[i - 1];

                double a1 = ((f2[i - 1] / (6.0f * h)));

                double a2 = ((f2[i] / (6.0f * h)));

                double b = ((y[i - 1] / h) - ((f2[i - 1] * h) / 6.0f));

                double c = ((y[i] / h) - ((f2[i] * h) / 6.0f));

                vector <double> coeficientes;

                coeficientes.push_back(a1);
                coeficientes.push_back(a2);
                coeficientes.push_back(b);
                coeficientes.push_back(c);

                int j = 1;

                cout << "\nEcuacion del Trazador Cubico en el intervalo [" << x[i - 1] << ", " << x[i] << "]" << endl
                     << "f" << j << "(x) = ";
              
                if (a1 != 0.0f){
                    cout << ((a1 < 0) ? " - " : "") << fabs(a1) << "*(" << x[i] << " - x)^3";
                }
                if (a2 != 0.0f){
                    cout << ((a2 < 0) ? " - " : " + ") << fabs(a2) << "*(x - " << x[i - 1] << ")^3";
                }

                cout << ((b < 0) ? " - " : " + ") << fabs(b) << "*(" << x[i] << " - x)";
                cout << ((c < 0) ? " - " : " + ") << fabs(c) << "*(x - " << x[i - 1] << ")" << endl;

                j++;

                return coeficientes;

            }

            /**
             * @brief Interpolar los coeficientes de cada subintervalo segun el intervalo de x_inicial a x_final
             * @param x_inicial Punto inicial del intervalo
             * @param x_final Punto final del intervalo
            */
            void trazadores(double x_inicial, double x_final){

                int n = x.size(); /*!< Numero de datos*/
                vector <double> coeficientes;
                int intervalos = n - 1; /*!< Numero de intervalos*/

                cout << "\nSegundas Derivadas: " << endl;
                for (size_t i = 0; i < f2.size(); i++){
                    cout << "f''(" << x[i] << ") = " << setprecision(6) << f2[i] << endl;
                }

                // Encontrar el primer intervalo en donde se encuentra x_inicial
                int i = 1;
                while (i < intervalos && x[i] < x_inicial) {
                    i++;
                }

                // Iterar sobre los intervalos hasta que se llegue al final del intervalo
                while (i <= intervalos && x[i - 1] <= x_final) { 
                    coeficientes = interpolar_trazador(x_inicial);
                    x_inicial = x[i++];
                }
                
            }
                            
            private:
                vector <double> x;
            vector <double> y;
            vector <double> f2;
            vector <double> calcular_f2(){
                
                vector <double> resultado;

                size_t n = x.size();
                size_t intervalos = n - 1;
                size_t i;

                vector <vector <double>> m(intervalos - 1);
                for (i = 0; i < intervalos - 1; i++){
                    m[i].resize(n);
                }

                for ( i = 1; i < intervalos; i++ ){

                    size_t fila = i - 1;

                    // * Primer coeficiente
                    if (i > 1) {
                        // Los puntos interiores tienen f''(xi-1)
                        m[fila][i - 1] = (x[i] - x[i - 1]);
                    };
                    
                    // * Segundo coeficiente
                    m[fila][i] = 2.0f * (x[i + 1] - x[i - 1]);

                    // * Tercer coeficiente
                    if (i < intervalos - 1){
                        // Los puntos interiores tienen f''(xi+1)
                        m[fila][i + 1] = (x[i + 1] - x[i]);
                    };
                    
                    double ci_1 = (6/(x[i + 1] - x[i])) * (y[i + 1] - y[i]);
                    double ci_2 = (6/(x[i] - x[i - 1])) * (y[i - 1] - y[i]);
                    double ci = ci_1 + ci_2;
                    m[fila][intervalos] = ci;
                    
                };

                // Eliminar los 0 al comienzo de las filas de la matriz
                for (i = 0; i < intervalos - 1; i++){
                    m[i].erase(m[i].begin());
                }

                resultado = gauss(m);

                vector <double> c;

                c.push_back(0);

                for (i = 0; i < resultado.size(); i++){
                    c.push_back(resultado[i]);
                }

                c.push_back(0);

                return c;
        };
    };
}

#endif