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
            
            /**
             * @brief Interpola el valor de x_int utilizando todos los datos
             * @param x_int Valor de x a interpolar
             * @return Valor interpolado
            */
            double interpolar(double x_int){
                return interpolar(x_int, 0, x.size() - 1);
            }

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
                    }
                }
                
                // Mostrar los factores del polinomio
                for(i = 0; i < b.size(); i++){
                    std::cout << "b" << i << " = " << b[i] << std::endl;
                }

                return oss.str();
            };

            /**
             * @brief Interpola el valor de x_int con Newton utilizando un polinomio del grado especifico
             * @param x_int Valor de x a interpolar sobre el cual se calcula el polinomio p(x)
             * @param pos_inicial Posicion inicial del intervalo
             * @param pos_final Posicion final del intervalo
             * @return Valor interpolado
            */
            double interpolar(double x_int, int pos_inicial, int pos_final){

                int n = x.size();

                // Validar que los coeficientes existan, y que pos_inicial y pos_final esten dentro del rango
                if (b.size() == 0 || pos_inicial < 0 || pos_final >= n || pos_inicial > pos_final) {return NAN;}

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
             * @brief Interpola el valor de x_int utilizando un polinomio del grado especifico
             * @param x_int Valor de x a interpolar sobre el cual se calcula el polinomio p(x)
             * @param grado Grado del polinomio p(x)
             * @return Valor interpolado
            */
            double interpolar(double x_int, int grado){              

                // 0. Encontrar las posiciones de x_int dentro de los datos
                int pos = lower_bound(x.begin(), x.end(), x_int) - x.begin();
                int n = x.size() - 1;
                            
                if (x_int < x[0] || x_int > x[n]) {
                    // x_int está por fuera del rango de los datos
                    return NAN;
                }

                // 1. Si n_puntos es par (tomar un solo intervalo) 
                if (grado % 2 == 0) {
                    int pos_Inicial = max(0, pos - grado / 2);
                    int pos_Final = min(n, pos_Inicial + grado);
                    return interpolar(x_int, pos_Inicial, pos_Final);
                }
                            
                // 2. Si n_puntos es impar (tomar dos intervalos)
                int pos_Inicial = max(0, pos - grado / 2);
                int pos_Final = min(n, pos_Inicial + grado);
                            
                double yInt_1 = interpolar(x_int, pos_Inicial, pos_Final);
                double error_int_1 = abs(calcular_error_interpolacion(x_int, pos_Inicial, pos_Final));
                            
                if (pos_Final + 1 <= n) {
                    double yInt_2 = interpolar(x_int, pos_Inicial + 1, pos_Final + 1);
                    double error_int_2 = abs(calcular_error_interpolacion(x_int, pos_Inicial + 1, pos_Final + 1));

                    if (abs(error_int_1) < abs(error_int_2)) {
                        return yInt_1;
                    }
                    return yInt_2;
                }
                
                return yInt_1;
            }

            /**
             * @brief Calcula el error de interpolacion usando pos_inicial y pos_final
             * @param x_int Valor de x a interpolar
             * @param pos_inicial Posicion inicial del intervalo
             * @param pos_final Posicion final del intervalo
             * @return Error de interpolacion
            */
            double calcular_error_interpolacion(double x_int, int pos_inicial, int pos_final){

                double valor_interpolado = interpolar(x_int, pos_inicial, pos_final);

                int pos = lower_bound(x.begin(), x.end(), x_int) - x.begin();
                double valor_real = y[pos];

                return valor_interpolado - valor_real;
            }

            /**
             * @brief Calcula el error de interpolacion el grado especificado
             * @param x_int Valor de x a interpolar
             * @param grado Grado del polinomio p(x)
             * @return Error de interpolacion 
            */
            double calcular_error_interpolacion(double x_int, int grado){

                double valor_interpolado = interpolar(x_int, grado);

                int pos = lower_bound(x.begin(), x.end(), x_int) - x.begin();
                double valor_real = y[pos];

                return valor_interpolado - valor_real;
                
            }

            /**
             * @brief Calcula el error de interpolacion usando todos los datos
             * @param x_int Valor de x a interpolar
             * @return Error de interpolacion
            */
            double calcular_error_interpolacion(double x_int){
                    
                double valor_interpolado = interpolar(x_int);

                int pos = lower_bound(x.begin(), x.end(), x_int) - x.begin();
                double valor_real = y[pos];

                return valor_interpolado - valor_real;
            }

            

            /**
             * @brief Metodo estatico para calcular los coeficientes del polinomio
             * @param x Variable independiente
             * @param y Variable dependiente
             * @param b Vector de coeficientes
             * @return Vector de coeficientes
            */
            vector<double> static calcular_coeficientes(vector<double> x, vector<double> y){
                
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
                return f[0];

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