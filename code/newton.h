/**
 * @file
 * @brief Interpolación mediante el método de diferencias divididas de Newton
 * @author Carlos Mario Perdomo Ramos <cmperdomo@unicauca.edu.co>
 * @author Daniel Fernando Solarte Ortega <dfsolarte@unicauca.edu.co>
*/

#ifndef NEWTON_H
#define NEWTON_H

#include <vector>
#include <cmath>
#include <sstream>
#include <string>

using std::vector;
using std::string;
using std::ostringstream;

namespace interpolacion{

    /**
     * @brief Metodo de Diferencias Divididas de Newton
    */
    class newton{

        public:
            /**
             * @brief Crea una instancia de Newton
             * @param p_x Variable independiente
             * @param p_y Variable dependiente
            */
            newton(vector<double> p_x, vector<double> p_y):x(p_x),y(p_y){
                calcular_coeficientes();
            }
            
            // TODO: Implementar metodo publico static para calcular los coeficientes y usarlo en spline

            /**
             * @brief Construye y retorna el polinomio interpolante
            */
            string polinomio(){
                ostringstream oss; // Flujo de salida de string

                size_t i, j;

                oss << b[0];
                for(i = 1; i < b.size(); i++){
                    oss << ((b[i] < 0)? " - ":" + ")
                        << fabs(b[i]) << " ";
                    for (j = 0; j < i; j++){
                        oss << "(x - " << x[j] << ") ";
                    };
                }
                
                // Mostrar los factores del polinomio
                for(i = 0; i < b.size(); i++){
                    cout << "b" << i << " = " << b[i] << endl;
                }

                return oss.str();
            };

            /**
             * @brief Calcula el valor de y interpolado
             * @param x_int Valor de x
             * @return y interpolado
            */
            double interpolar(double x_int){
                if (b.size() == 0) {return NAN;}
                double f = b[0];

                size_t i, j;

                for(i = 1; i < b.size(); i++){
                    double prod = 1.0f;
                    for (j = 0; j < i; j++){
                        prod *= (x_int - x[j]);
                    };

                    f += b[i] * prod;
                }

                return f;
            }

            /**
             * @brief Metodo estatico para calcular los coeficientes del polinomio
             * @param x Variable independiente
             * @param y Variable dependiente
             * @param b Vector de coeficientes
            */
            void static calcular_coeficientes(vector<double> x, vector<double> y, vector<double> &b){
                
                size_t i,j;
                size_t n = x.size();
                
                // Matriz: vector de vectores
                vector <vector <double>> f(n);

                for (i = 0; i<n; i++){
                    f[i].resize(n - i);
                }

                // Llenar la primera columna
                for (i = 0; i < n; i++){
                    f[i][0] = y[i];
                }

                // Llenar la matriz con las diferencias divididas
                for (j = 1; j < n; j++){
                    for (i = 0; i < n - j; i++){
                        f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i]);
                    }
                }

                // Tomar los coeficientes de la primera fila de la matriz
                b = f[0];

            } 

        private:
            /** @brief Calcula los coeficientes del polinomio */
            void calcular_coeficientes(){
                
                size_t i,j;
                size_t n = x.size();
                
                // Matriz: vector de vectores
                vector <vector <double>> f(n);

                for (i = 0; i<n; i++){
                    f[i].resize(n - i);
                }

                // Llenar la primera columna
                for (i = 0; i < n; i++){
                    f[i][0] = y[i];
                }

                // Llenar la matriz con las diferencias divididas
                for (j = 1; j < n; j++){
                    for (i = 0; i < n - j; i++){
                        f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i]);
                    }
                }

                // Tomar los coeficientes de la primera fila de la matriz
                b = f[0];

            }
            vector <double> x; /*!< Variable independiente */
            vector <double> y; /*!< Variable dependiente */
            vector <double> b; /*!< Coeficientes b0, b1, ... del polinomio */

    };
}

#endif