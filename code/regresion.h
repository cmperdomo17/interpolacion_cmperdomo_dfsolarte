#ifndef REGRESION_H
#define REGRESION_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <cmath>

#include "util.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;
using std::setw;
using std::setprecision;
using std::left;
using std::right;
using std::to_string;

using util::gauss;

namespace regresion{

    /**
     * @brief Repite una cadena de caracteres
     * @param s Cadena a repetir
     * @param n Cantidad de repeteciones
     * @return Concatenacion de s n veces
    */
    string str_repeat(string s, int n){
        string ret = "";
        for(int i = 0; i < n; i++){
            ret += s;
        }
        return ret;
    }

    int getCenteredWidth(size_t label_length, int value) {
        size_t value_length = to_string(value).length();
        return (label_length - value_length) / 2 + value_length;
    }



    /**
     * @brief Solucion mediante Regresion Lineal Simple
    */
    struct solucion_lineal{
        double b0 = NAN; /*!<Termino independiente de la recta */
        double b1 = NAN; /*!< Coeficiente de X */
        double st; /*!< Sumatoria de la diferencia cuadratica entre el valor medido y el promedio */
        double sy; /*!< Desviacion estandar */
        double sr; /*!< Sumatoria de la diferencia cuadratica entre cada y con el y estimado */
        double syx; /*!< Error estandar de aproximacion */
        double r2; /*!< Coeficiente de determinacion */
        size_t n; /*!< Numero de datos */

        /**
         *
         */
         void imprimir(){

            string aceptable = (syx < sy) ? "La aproximacion se considera aceptable" : "La aproximacion NO se considera aceptable";

             cout << "Recta de regresion: \n"
                 << "\ny = " << b1 << " * x "
                 << ((b0 >= 0.0f)? " + " : " - ")
                 << fabs(b0)
                 << "\n"
                 << endl   
                 << "Desviacion estandar: "
                 << sy
                 << "\n"
                 << endl
                 << "Error estandar de aproximacion: "
                 << syx
                 << "\n"
                 << aceptable
                 << endl
                 << "Coeficiente de determinacion: "
                 << r2
                 << endl;
         }
    };

    /**
     * @brief Solución Potencia
    */
    struct solucion_potencia
    {
        double c = 0.0f;        /*!<Coeficiente de la potencia*/
        double a = 0.0f;        /*!<Factor del exponente de la potencia*/
        solucion_lineal lineal; /*!<Regresion de los datos linealizados*/
        /**
         * @brief Impresion de la regresion potencia
        */
        void imprimir()
        {
            string aceptable = (lineal.syx < lineal.sy) ? "La aproximacion se considera aceptable" : "La aproximacion NO se considera aceptable";
            cout << "Recta de regresion" << endl
                 << "y= " << c << " *x^" << a
                 << endl
                 << "Desviacion estandar: " 
                 << lineal.sy 
                 << endl
                 << "Error estandar de aproximacion: " 
                 << lineal.syx
                 << endl
                 << aceptable
                 << endl
                 << "Coeficiente de determinacion: "  
                 << lineal.r2 
                 << endl;
        }
    };

    /**
     * @brief Solución exponente
    */
    struct solucion_exponencial
    {
        double c = 0.0f;        /*!<Coeficiente de la potencia*/
        double a = 0.0f;        /*!<Factor del exponente de la potencia*/
        solucion_lineal lineal; /*!<Regresion de los datos linealizados*/
        /**
         * @brief Impresion de la regresion potencia
        */
        void imprimir()
        {
            string aceptable = (lineal.syx < lineal.sy) ? "La aproximacion se considera aceptable" : "La aproximacion NO se considera aceptable";
            cout << "Recta de regresion" << endl
                 << "y= " << c << " *e^" << a << "*x"
                 << endl
                 << "Desviacion estandar: " 
                 << lineal.sy 
                 << endl
                 << "Error estandar de aproximacion: " 
                 << lineal.syx
                 << endl
                 << aceptable
                 << endl
                 << "Coeficiente de determinacion: "  
                 << lineal.r2 
                 << endl;
        }
    };

    /**
     * @brief Solución cuadrática
    */
    struct solucion_cuadratica{

        double a0 = 0.0f; /*!< Termino independiente del polinomio cuadratico */
        double a1 = 0.0f; /*!< Coeficiente de x del polinomio cuadratico */
        double a2 = 0.0f; /*!< Coeficiente de x^2 del polinomio cuadratico */
        double st = 0.0f; /*!< Sumatoria de la diferencia cuadratica entre el valor medido y el promedio */
        double sr = 0.0f; /*!< Sumatoria de la diferencia cuadratica entre cada y con el y estimado */
        double sy = 0.0f; /*!< Desviacion estandar */
        double syx = 0.0f; /*!< Error estandar de aproximacion */
        double r2 = 0.0f; /*!< Coeficiente de determinacion */
        size_t n; /*!< Numero de datos */

        /**
         * @brief Impresion de la regresion cuadratica
        */
        void imprimir(){

            string aceptable = (syx < sy) ? "La aproximacion se considera aceptable" : "La aproximacion NO se considera aceptable";
            
            cout << "Polinomio de regresion: \n"
            << "y = " << a2 << " * x^2 "
            << ((a1 >= 0.0f)? " + " : " - ")
            << fabs(a1) << "x"
            << ((a0 >= 0.0f)? " + " : " - ")
            << fabs(a0)
            << "\n"
            << endl
            << "Desviacion estandar: "
            << sy
            << "\n"
            << endl
            << "Error estandar de aproximacion: "
            << syx
            << "\n"
            << endl
            << "Coeficiente de determinacion: "
            << r2
            << endl;
        }
    };

    /**
     * Imprime una tabla de datos
     * @param X  variable independiente
     * @param y variable dependiente
     * @param x_label Etiqueta de la variable independiente
     * @param y_label Etiqueta de la variable dependiente
     */
    void imprimir_tabla(vector<double> x, vector<double> y, string x_label = "", string y_label = "") {

        if (x_label.empty()) {
            x_label = "X";
        }

        if (y_label.empty()) {
            y_label = "Y";
        }

        size_t x_width = x_label.length() + 4;
        size_t y_width = y_label.length() + 4;

        cout << "\n";
        cout << str_repeat("=", x_width + y_width + 2) << endl;

        cout << "  " << setw(x_width) << left << x_label << setw(y_width/4) << right << y_label << endl;

        cout << str_repeat("=", x_width + y_width + 2) << endl;

        for (size_t i = 0; i < x.size(); i++) {
            int x_space = getCenteredWidth(x_width, x[i]);
            int y_space = getCenteredWidth(y_width, y[i]);

            cout << setw(x_space) << right << x[i] << setw(x_width - x_space + y_space) << right << y[i] << endl;
        }
        cout << str_repeat("=", x_width + y_width + 2) << endl;
        cout << "\n";

        }

    class lineal_simple{
    public:
        /**
         *
         */
        lineal_simple(vector<double> p_x, vector<double> p_y):x(p_x), y(p_y){
        }

        /**
         * @brief Calcula la recta de regresion lineal
         * @return  Recta de regresion lineal
         */
        solucion_lineal calcular(){
            solucion_lineal sol;

            double sum_xy = 0.0f, sum_x = 0.0f, sum_x2 = 0.0f, sum_y = 0.0f,
                x_prom , y_prom;

            sol.n = x.size();

            for(size_t i = 0; i < sol.n; i++){
                sum_xy += x[i] * y[i];
                sum_x += x[i];
                sum_y += y[i];
                sum_x2 += x[i] * x[i];
            }

            if(sol.n > 0){
                x_prom = sum_x / sol.n;
                y_prom = sum_y / sol.n;

                sol.st = 0.0f;
                for(size_t i = 0; i < sol.n; i++){
                    sol.st += pow(y[i] - y_prom, 2.0f);
                }

                if (sol.n > 1){
                    sol.sy = sqrt(sol.st / (double)(sol.n - 1));
                }
            } else {
                x_prom = NAN;
                y_prom = NAN;
            }

            //Calcular b1 y b0
            sol.b1 = (sum_xy - (y_prom * sum_x)) / (sum_x2 - (x_prom * sum_x));
            sol.b0 = y_prom - (sol.b1 * x_prom);

            sol.sr = 0.0f;
            for(size_t i = 0; i < sol.n; i++){
                sol.sr += pow(y[i] - ((sol.b1 * x[i]) + sol.b0), 2.0f);
            }

            if(sol.n > 2){
                sol.syx = sqrt(sol.sr / (double)(sol.n - 2));
            }

            sol.r2 = (sol.st - sol.sr) / sol.st;
            return sol;
        }

    private:
        vector<double> x; /*!< Variable independiente */
        vector<double> y; /*!< Variable dependiente */
    };


    /**
     * @brief Regresion linealizada mediante la funcion potencia
    */
    class potencia{
    
    public:
        /**
         * @brief Crea una nueva instancia de regresion mediante funcion potencia
        */
        potencia(vector<double> p_x, vector<double> p_y):x(p_x), y(p_y){
        }

        /**
         * @brief Calcula la regresion linealizada mediante la funcion potencia
         * @return Solucion linealizada
        */
        solucion_potencia calcular(){

            solucion_potencia sol;

            vector<double> X(x);
            vector<double> Y(y);

            for(unsigned int i = 0; i < X.size(); i++){
                X[i] = log10(X[i]);
                Y[i] = log10(Y[i]);
            }

            // Crear un modelo de regresion lineal con los datos transformados
            lineal_simple ls(X, Y);

            // Calcular la regresion lineal con los datos transformados
            sol.lineal = ls.calcular();

            // Calcular A
            sol.a = sol.lineal.b1;

            // Calcular C
            sol.c = pow(10.0f, sol.lineal.b0);

            return sol;

        }

    private:
        vector<double> x; /*!< Variable independiente */
        vector<double> y; /*!< Variable dependiente */   
    };

    /**
     * @brief Regresion linealizada mediante la funcion potencia
    */
    class exponencial{
    
    public:
        /**
         * @brief Crea una nueva instancia de regresion mediante funcion exponencial
        */
        exponencial(vector<double> p_x, vector<double> p_y):x(p_x), y(p_y){
        }

        /**
         * @brief Calcula la regresion linealizada mediante la funcion potencia
         * @return Solucion linealizada
        */
        solucion_exponencial calcular(){

            solucion_exponencial sol;

            vector<double> X(x);
            vector<double> Y(y);

            for(unsigned int i = 0; i < Y.size(); i++){
                Y[i] = log(Y[i]);
            }

            // Crear un modelo de regresion lineal con los datos transformados
            lineal_simple ls(x, Y);

            // Calcular la regresion lineal con los datos transformados
            sol.lineal = ls.calcular();

            // Calcular A
            sol.a = sol.lineal.b1;

            // Calcular C
            sol.c = exp(sol.lineal.b0);

            return sol;

        }

    private:
        vector<double> x; /*!< Variable independiente */
        vector<double> y; /*!< Variable dependiente */   
    };

    /**
     * @brief Regresion cuadratica
    */
    class cuadratica{

        public:
            cuadratica(vector<double> p_x, vector<double> p_y):x(p_x), y(p_y){
            };
            /**
             * @brief Calcula el polinomio de regresion de grado 2
             * @return Polinomio de solucion.
            */
            solucion_cuadratica calcular(){
                
                solucion_cuadratica sol;
                size_t i;

                double sum_x{0}, sum_x2{0}, sum_x3{0}, sum_x4{0};
                double sum_y{0}, sum_xy{0}, sum_x2y{0};
                double y_prom;

                sol.n = x.size();

                // La regresion cuadratica se calcula con al menos 4 puntos
                if (sol.n <= 3){
                    return sol;
                }

                for(i=0; i<sol.n; i++){
                    sum_x += x[i];
                    double x2 = pow(x[i], 2.0f);
                    sum_x2 += x2;
                    sum_x3 += pow(x[i], 3.0f);
                    sum_x4 += pow(x[i], 4.0f);
                    sum_y += y[i];
                    sum_xy += x[i] * y[i];
                    sum_x2y += x2 * y[i];
                }

                // Calcular y_prom
                y_prom = sum_y / (double)sol.n;

                

                /*
                cout << "sum_x: " << sum_x << endl;
                cout << "sum_x2: " << sum_x2 << endl;
                cout << "sum_x3: " << sum_x3 << endl;
                cout << "sum_x4: " << sum_x4 << endl;
                cout << "sum_y: " << sum_y << endl;
                cout << "sum_xy: " << sum_xy << endl;
                cout << "sum_x2y: " << sum_x2y << endl;
                */

                vector <vector<double>> m {
                    {(double)sol.n, sum_x, sum_x2, sum_y},
                    {sum_x, sum_x2, sum_x3, sum_xy},
                    {sum_x2, sum_x3, sum_x4, sum_x2y}
                };
                
                // Hallar a0, a1 y a2 mediante eliminacion de Gauss
                vector <double> coef = gauss(m);

                //Imprimir los coeficientes

                sol.a0 = coef[0];
                sol.a1 = coef[1];
                sol.a2 = coef[2];

                // Calcular st
                sol.st = 0.0f;
                for(i=0; i<sol.n; i++){
                    sol.st += pow(y[i] - y_prom, 2.0f);
                }

                // Calcular sy
                sol.sy = sqrt(sol.st / (double)(sol.n - 1));  

               	//Calcular sr
			    sol.sr = 0.0f;
                for(size_t i = 0; i<sol.n; i++) {
                    sol.sr += pow(y[i] - sol.a0 - (sol.a1*x[i]) - (sol.a2*pow(x[i], 2.0f)) , 2.0f);
                }
			
                //Calcular sxy
                if(sol.n > 2) {
                    sol.syx = sqrt(sol.sr/(double)(sol.n - 3));
                }
			
                //Calcular r2
                sol.r2 = (sol.st - sol.sr)/sol.st;
                
                return sol;

            }
        private:
            vector<double> x; /*!< Variable independiente */
            vector<double> y; /*!< Variable dependiente */
    };
}
#endif