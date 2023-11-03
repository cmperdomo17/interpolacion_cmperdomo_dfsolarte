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
#include <iostream>

using std::vector;
using std::string;
using std::ostringstream;
using std::isnan;
using std::endl;
using std::cout;
using std::lower_bound;
using std::fabs;
using std::invalid_argument;

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
                oss << "\n\nFactores del polinomio: " << endl;
                oss << endl;
                for(i = 0; i < b.size(); i++){          
                    oss << "b" << i << " = " << b[i] << endl;
                }

                return oss.str();
            };

            /**
             * @brief Interpola el valor de x_int con Newton utilizando un polinomio del grado especifico
             * @param x_int Valor de x a interpolar sobre el cual se calcula el polinomio p(x)
             * @param pos_Inicial Posicion inicial del intervalo
             * @param pos_Final Posicion final del intervalo
             * @return Valor interpolado
            */
            double interpolar(double x_int, int pos_Inicial, int pos_Final){

                int n = x.size();

                // Validar que los coeficientes existan, y que pos_Inicial y pos_Final esten dentro del rango
                if (b.size() == 0 || pos_Inicial < 0 || pos_Final >= n || pos_Inicial > pos_Final) {return NAN;}

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
                //Validar que x_int este dentro del rango de x
                if (x_int < x[0] || x_int >= x[x.size() - 1]){
                    return NAN;
                }

                int n_puntos = grado + 1,
                        n = x.size(),
                        pos_ant = 0,
                        pos_sig = x.size() - 1,
                        pos_Inicial = 0,
                        pos_Final = n - 1,
                        pos_Inicial_aux,
                        pos_Final_aux;
                //Validar que el grado no sea mayor a n_puntos
                if (grado < 0 || grado >= n_puntos){
                    throw invalid_argument("Grado invalido");
                }

                //Encontrar la posicion anterior y siguiente de x_int dentro de los datos originales
                for (size_t i = 0; i < x.size(); i++) {
                    if (x[i] <= x_int) {
                        pos_ant = i;
                    } else if (x[i] > x_int) {
                        pos_sig = i;
                        break;
                    }
                }

                pos_Inicial = pos_ant - grado / 2;
                pos_Final = pos_sig + grado / 2;
                //Si n_puntos es par
                if (n_puntos % 2 == 0) {
                    return interpolar(x_int, pos_Inicial, pos_Final);
                } else { 
                    // Si n_puntos es impar
                    pos_Inicial_aux = pos_Inicial + 1;
                    pos_Final_aux = pos_Final - 1;

                    // Validar Posiciones
                    if (pos_Inicial_aux >= pos_Final_aux || pos_Final_aux >= n || pos_Inicial_aux >= n || pos_Inicial < 0 || pos_Final >= n) {
                        return interpolar(x_int, pos_Inicial, pos_Final);
                    }

                    double y_int_1 = interpolar(x_int, pos_Inicial, pos_Final);
                    double y_int_2 = interpolar(x_int, pos_Inicial_aux, pos_Final_aux);

                    if (isnan(y_int_1)) {
                        return y_int_2;
                    } else if(isnan(y_int_2)) {
                        return y_int_1;
                    }

                    // y_int_1 o y_int_2 son diferente de nan
                    //Sacar los datos de x en el intervalo pos_Inicial_aux, pos_Final_aux con el dato adicional
                    // con un for del x grande (pos_Inicial) sacr los datos a x1 (0), x1 es un subvector de x que tiene desde x[pos_Inicial_aux] hasta x[pos_Final_aux]
                    vector<double> x1 (x.begin() + pos_Inicial, x.begin() + pos_Final);

                    //Sacar los datos de y en el intervalo pos_Inicial_aux, pos_Final_aux con el dato adicional
                    vector<double> y1 (y.begin() + pos_Inicial, y.begin() + pos_Final);

                    vector<double> F1 = calcular_coeficientes(x1, y1);

                    //Calcular el error
                    double prod_1 = 1.0f; //Ultimo coeficiente de F1

                    //Quitar el dato adicional del fin
                    F1.erase(F1.end());

                    //Calcular la productoria de R * (x_int - x1[0]) * (x_int - x1[1]) * ... * (x_int - x1[n_puntos - 1]) sin tener en cuenta el dato adicional
                    for (size_t i = 0; i < F1.size(); i++) {
                        prod_1 *= (x_int - x1[i]);
                    }
                    double error_int_1 =  F1[F1.size() - 1] * prod_1 ;

                    //Imprimir primer intervalo
                    //Sacar los datos de x en el intervalo pos_Inicial_aux, pos_Final_aux con el dato adicional
                    cout << "\nPrimer intervalo: " << endl;
                    cout << "   Posicion Inicial: " << pos_Inicial << ", Posicion Final: " << pos_Final << endl;
                    cout << "   Error 1 (R1): " << error_int_1 << endl;

                    vector<double> x2 (x.begin() + pos_Inicial_aux, x.begin() + pos_Final_aux);

                    //Sacar los datos de y en el intervalo pos_Inicial_aux, pos_Final_aux con el dato adicional
                    vector<double> y2 (y.begin() + pos_Inicial_aux, y.begin() + pos_Final_aux);

                    vector<double> F2 = calcular_coeficientes(x2, y2);

                    //Calcular el error
                    double prod_2 = 1.0f; //Ultimo coeficiente de F1

                    //Quitar el dato adicional del inicio
                    F2.erase(F2.begin());

                    //Calcular la productoria de R * (x_int - x1[0]) * (x_int - x1[1]) * ... * (x_int - x1[n_puntos - 1]) sin tener en cuenta el dato adicional
                    for (size_t i = 0; i < F2.size(); i++) {
                        prod_2 *= (x_int - x1[i]);
                    }
                    double error_int_2 =  F2[F2.size() - 1] * prod_2;

                    //Imprimir segundo intervalo
                    //Sacar los datos de x en el intervalo pos_Inicial_aux, pos_Final_aux con el dato adicional
                    cout << "\nSegundo intervalo: " << endl;
                    cout << "   Posicion Inicial: " << pos_Inicial_aux << ", Posicion Final: " << pos_Final_aux << endl;
                    cout << "   Error 2 (R2): " << error_int_2 << endl;

                    //Imprimir el error superior e inferior de los intervalos

                    if (fabs(error_int_1) < fabs(error_int_2)) {
                        return y_int_2;
                    } else {
                        return y_int_1;
                    }
                }
                return NAN;
            }

            /**
             * @brief Calcula el error de interpolacion usando pos_Inicial y pos_Final
             * @param x_int Valor de x a interpolar
             * @param pos_Inicial Posicion inicial del intervalo
             * @param pos_Final Posicion final del intervalo
             * @return Error de interpolacion
            */
            double calcular_error_interpolacion(double x_int, int pos_Inicial, int pos_Final){

                double valor_interpolado = interpolar(x_int, pos_Inicial, pos_Final);

                int pos = lower_bound(x.begin(), x.end(), x_int) - x.begin();
                double valor_real = y[pos];

                return valor_interpolado - valor_real;
            }

            /**
             * @brief Calcula el error de interpolacion usando el grado especificado
             * @param x_int Valor de x a interpolar
             * @param grado Grado del polinomio p(x)
             * @return Error de interpolacion 
            */
            double calcular_error_interpolacion(double x_int, double valor_interpolado){

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