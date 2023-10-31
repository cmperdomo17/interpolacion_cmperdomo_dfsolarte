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
                // Recorrer los intervalos hasta encontrar el intervalo i
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

                // TODO 2: Calcular F2
                // * Calcular F2: F2 = gauss(M)
                // * Insertar 0 al inicio y al final de F2
                // * (f2 en los extremos vale 0)

                // matriz con los ceros al comienzo
                for (i=0; i <m.size(); i++){
                    for (size_t j=0; j<m[i].size(); j++){
                        cout << m[i][j] << " ";
                    }
                    cout << endl;
                }

                // quitar los 0 al comienzo de las filas de la matriz
                for (i = 0; i < intervalos - 1; i++){
                    m[i].erase(m[i].begin());
                }

                cout << "\n" << endl;

                // matriz sin los 0 al comienzo de las filas
                for (i=0; i <m.size(); i++){
                    for (size_t j=0; j<m[i].size(); j++){
                        cout << m[i][j] << " ";
                    }
                    cout << endl;
                }
    
                resultado = gauss(m);

                for (i=0; i<resultado.size(); i++){
                   cout << resultado[i] << endl;
                }

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