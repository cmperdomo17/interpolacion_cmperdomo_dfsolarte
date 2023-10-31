#include <iostream>
#include <vector>
#include "lagrange.h"
#include "newton.h"
#include "util.h"
#include "spline3.h"
#include "regresion.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;

using interpolacion::newton;
using interpolacion::lagrange;
using util::imprimir_tabla;
using interpolacion::spline3;

using regresion::lineal_simple;
using regresion::exponencial;
using regresion::potencia;
using regresion::cuadratica;

using regresion::solucion_lineal;
using regresion::solucion_potencia;
using regresion::solucion_exponencial;
using regresion::solucion_cuadratica;

/** 
 * @brief Parte 1: Regresion
*/
void caso_regresion();

/** 
 * @brief Parte 2.a: Interpolacion Diferencias Divididas de Newton
*/
void caso_interpolacion_newton();

/** 
 * @brief Parte 2.b: Interpolacion Trazadores Cubicos
*/
void caso_interpolacion_spline3();

/** 
 * @brief Parte 2.c: Interpolacion Lagrange
*/
void caso_interpolacion_lagrange();

int main() {
    int opcion;

    do {
        cout << "\n ------ MENU PRINCIPAL ------\n" << endl;
        cout << "1. Parte 1: Regresion " << endl;
        cout << "2. Parte 2: Interpolacion" << endl;
        cout << "0. Salir" << endl;
        cout << "\nIngrese una opcion: ";
        cin >> opcion;

        switch (opcion) {
            case 1:
                caso_regresion();
                break;
            case 2:
                do {
                    cout << "\n ------ MENU INTERPOLACION ------\n" << endl;
                    cout << "1. Caso 1: Interpolacion por diferencias divididas de Newton" << endl;
                    cout << "2. Caso 2: Interpolacion por Trazadores Cubicos" << endl;
                    cout << "3. Caso 3: Interpolacion por Lagrange" << endl;
                    cout << "4. Volver al menu principal" << endl;
                    cout << "0. Salir" << endl;
                    cout << "\nIngrese una opcion: ";
                    cin >> opcion;

                    switch (opcion) {
                        case 1:
                            caso_interpolacion_newton();
                            break;
                        case 2:
                            caso_interpolacion_spline3();
                            break;
                        case 3:
                            caso_interpolacion_lagrange();
                            break;
                        case 4:
                            cout << "Volviendo al menu principal..." << endl;
                            break;
                        case 0:
                            cout << "Saliendo..." << endl;
                            break;
                        default:
                            cout << "Opcion invalida" << endl;
                            break;
                    }
                } while (opcion != 0 && opcion != 4);
                break;
            case 0:
                cout << "Saliendo..." << endl;
                break;
            default:
                cout << "Opcion invalida" << endl;
                break;
        }
    } while (opcion != 0);

    return 0;
}

// Casos de regresion

void caso_regresion_lineal(vector <double> x, vector <double> y, string title, string x_label, string y_label){
    // Imprime el titulo del caso
    cout << title << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, x_label, y_label);
    lineal_simple ls(x,y);
    solucion_lineal sol = ls.calcular();

    sol.imprimir();
}

void caso_regresion_potencia(vector <double> x, vector <double> y, string title, string x_label, string y_label){
    // Imprime el titulo del caso
    cout << title << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, x_label, y_label);
    potencia reg_pot(x,y);
    solucion_potencia sol = reg_pot.calcular();

    sol.imprimir();
}

void caso_regresion_exponencial(vector <double> x, vector <double> y, string title, string x_label, string y_label){
    // Imprime el titulo del caso
    cout << title << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, x_label, y_label);
    exponencial reg_exp(x,y);
    solucion_exponencial sol = reg_exp.calcular();

    sol.imprimir();
}

void caso_regresion_cuadratica(vector <double> x, vector <double> y, string title, string x_label, string y_label){
    // Imprime el titulo del caso
    cout << title << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, x_label, y_label);
    cuadratica reg_cua(x,y);
    solucion_cuadratica sol = reg_cua.calcular();

    sol.imprimir();
}

void caso_regresion(){

    vector <double> x = {4.0f, 5.0f, 8.0f, 11.0f, 14.0f, 17.0f, 19.0f, 23.0f, 25.0f, 26.0f, 29.0f, 35.0f, 48.0f};
    vector <double> y = {82.3f, 79.2f, 75.6f, 72.4f, 68.9f, 67.6, 66.8f, 65.1f, 64.8f, 63.5f, 62.6f, 61.4f, 60.3f};

    caso_regresion_lineal(x, y, "\nCaso 1: Regresion Lineal", "x", "y");
    caso_regresion_potencia(x, y, "\nCaso 2: Regresion Potencia", "x", "y");
    caso_regresion_exponencial(x, y, "\nCaso 3: Regresion Exponencial", "x", "y");
    caso_regresion_cuadratica(x, y, "\nCaso 4: Regresion Cuadratica", "x", "y");
}

// Casos de interpolacion

void caso_interpolacion_newton(){

    vector <double> x = {0.4f, 0.8f, 1.3f, 1.8f, 2.0f, 2.2f, 2.6f};
    vector <double> y = {1.452360f, 1.995632f, 2.719678f, 3.273019f, 3.359425f, 3.316678f, 2.669452f};

    // Instancia de Newton
    newton n(x, y);

    // Imprimir el polinomio
    cout << "\nPolinomio interpolante: \n\n" << n.polinomio() << endl;

    // Valor a interpolar
    double x_int;
    // Grado de interpolacion
    size_t grado;
    // Valor interpolado
    double y_int;
    // Error de interpolacion
    double error_int;
    // Error de interpolacion con las posiciones iniciales y finales
    double error_int_1; 

    cout << "Interpolacion por diferencias divididas de Newton" << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, "    X    ", "    Y    ");

    // Solicitar el valor a interpolar
    do{
        cout << "Ingrese el valor a interpolar: ";
        cin >> x_int;
    } while(x_int < x[0] || x_int > x[x.size() - 1]);

    // Solicitar el grado de interpolación
    do{
        cout << "Ingrese el grado del polinomio, <= : " << x.size() - 1 << ", 0 para usar todos los datos: ";
        cin >> grado;
    } while (grado > x.size());


    // Interpolar el valor ingresado por el usuario con todos los datos
    if (grado == 0){
        y_int = n.interpolar(x_int);
        // Calcula el error de interpolacion con todos los datos
        error_int = abs(n.calcular_error_interpolacion(x_int));
        // Mostrar la y interpolada
        cout << "\ny = " << setprecision(7) << y_int << endl;
        // Mostrar el error de interpolacion
        cout << "\nError de interpolacion: " << error_int << endl;
    } else {

        // Interpolar el valor ingresado por el usuario con el grado especificado
        y_int = n.interpolar(x_int, grado);
        // Calcula el error de interpolacion con el grado especificado
        error_int = abs(n.calcular_error_interpolacion(x_int, grado));

        // Mostrar el Primer Intervalo
        cout << "\nPrimer Intervalo: " << endl;
        cout << "\n - Posicion Inicial: " << 0 << ", Posicion Final: " << grado << endl;
        error_int_1 = abs(n.calcular_error_interpolacion(x_int, 0, grado));
        cout << " - R: " << error_int_1 << endl;

        // Mostrar el Segundo Intervalo
        cout << "\nSegundo Intervalo: " << endl;
        cout << "\n - Posicion Inicial: " << 1 << ", Posicion Final: " << grado + 1 << endl;
        error_int_1 = abs(n.calcular_error_interpolacion(x_int, 1, grado + 1));
        cout << " - R: " << error_int_1 << endl;

        // Calcular el índice del punto más cercano por encima de x_int
        size_t index_sup = 0;
        while (index_sup < x.size() && x[index_sup] <= x_int) {
            index_sup++;
        }

        // Calcular el índice del punto más cercano por debajo de x_int
        size_t index_inf = index_sup - 1;

        // Calcular el error superior (error en el punto por encima)
        double error_superior = abs(y_int - y[index_sup]);

        // Calcular el error inferior (error en el punto por debajo)
        double error_inferior = abs(y_int - y[index_inf]);

        double y_inferior = y[index_inf];
        double y_superior = y[index_sup];

        // Mostrar Y superior e inferior
        cout << "\nY (superior): " << setprecision(7) << y_superior << endl;
        cout << "Y (inferior): " << setprecision(7) << y_inferior << endl;
        cout << "Error (superior): " << setprecision(7) << error_superior << endl;
        cout << "Error (inferior): " << setprecision(7) << error_inferior << endl;

        // Mostrar la y interpolada
        cout << "\ny = " << setprecision(7) << y_int << endl;

        // Mostrar el mejor error de interpolacion
        if (abs(error_superior) < abs(error_inferior)){
            cout << "\nMejor error de interpolacion: " << error_superior << endl;
        } else {
            cout << "\nMejor error de interpolacion: " << error_inferior << endl;
        }
    }

}

void caso_interpolacion_spline3(){

    vector <double> x = {0.4f, 0.8f, 1.3f, 1.8f, 2.0f, 2.2f, 2.6f};
    vector <double> y = {1.452360f, 1.995632f, 2.719678f, 3.273019f, 3.359425f, 3.316678f, 2.669452f};

    // Instancia de Trazador Cubico
    spline3 s3(x, y);

    // Valor a interpolar
    double x_int;

    cout << "\nInterpolacion mediante el metodo de Trazadores Cubicos" << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, "    x    ", "   f(x)    ");
    // Solicitar el valor a interpolar
    do{
        cout << "Ingrese el valor a interpolar: ";
        cin >> x_int;
    } while(x_int < x[0] || x_int > x[x.size() - 1]);

    // Interpolar el valor ingresado por el usuario
    double y_int = s3.interpolar(x_int);

    cout << "\ny = " << setprecision(7) << y_int << endl;
    
}

void caso_interpolacion_lagrange(){

    vector <double> x = {0.4f, 0.8f, 1.3f, 1.8f, 2.0f, 2.2f, 2.6f};
    vector <double> y = {1.452360f, 1.995632f, 2.719678f, 3.273019f, 3.359425f, 3.316678f, 2.669452f};

    // Instancia de Lagrange
    lagrange l(x, y);

    // Imprimir el polinomio
    cout << "\nPolinomio interpolante: \n\n" << l.polinomio() << endl;

    // Valor a interpolar
    double x_int;
    // Grado de interpolacion
    size_t grado;
    // Valor interpolado
    double y_int;
    // Error de interpolacion
    double error_int;
    // Error de interpolacion con las posiciones iniciales y finales
    double error_int_1; 

    cout << "Interpolacion por diferencias divididas de Newton" << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, "    X    ", "    Y    ");

    // Solicitar el valor a interpolar
    do{
        cout << "Ingrese el valor a interpolar: ";
        cin >> x_int;
    } while(x_int < x[0] || x_int > x[x.size() - 1]);

    // Solicitar el grado de interpolación
    do{
        cout << "Ingrese el grado del polinomio, <= : " << x.size() - 1 << ", 0 para usar todos los datos: ";
        cin >> grado;
    } while (grado > x.size());


    // Interpolar el valor ingresado por el usuario con todos los datos
    if (grado == 0){
        y_int = l.interpolar(x_int);
        // Calcula el error de interpolacion con todos los datos
        error_int = abs(l.calcular_error_interpolacion(x_int));
        // Mostrar la y interpolada
        cout << "\ny = " << setprecision(7) << y_int << endl;
        // Mostrar el error de interpolacion
        cout << "\nError de interpolacion: " << error_int << endl;
    } else {

        // Interpolar el valor ingresado por el usuario con el grado especificado
        y_int = l.interpolar(x_int, grado);
        // Calcula el error de interpolacion con el grado especificado
        error_int = abs(l.calcular_error_interpolacion(x_int, grado));

        // Mostrar el Primer Intervalo
        cout << "\nPrimer Intervalo: " << endl;
        cout << "\n - Posicion Inicial: " << 0 << ", Posicion Final: " << grado << endl;
        error_int_1 = abs(l.calcular_error_interpolacion(x_int, 0, grado));
        cout << " - R: " << error_int_1 << endl;

        // Mostrar el Segundo Intervalo
        cout << "\nSegundo Intervalo: " << endl;
        cout << "\n - Posicion Inicial: " << 1 << ", Posicion Final: " << grado + 1 << endl;
        error_int_1 = abs(l.calcular_error_interpolacion(x_int, 1, grado + 1));
        cout << " - R: " << error_int_1 << endl;

        
        // Calcular el índice del punto más cercano por encima de x_int
        size_t index_sup = 0;
        while (index_sup < x.size() && x[index_sup] <= x_int) {
            index_sup++;
        }

        // Calcular el índice del punto más cercano por debajo de x_int
        size_t index_inf = index_sup - 1;

        // Calcular el error superior (error en el punto por encima)
        double error_superior = abs(y_int - y[index_sup]);

        // Calcular el error inferior (error en el punto por debajo)
        double error_inferior = abs(y_int - y[index_inf]);

        double y_inferior = y[index_inf];
        double y_superior = y[index_sup];

        // Mostrar Y superior e inferior
        cout << "\nY (superior): " << setprecision(7) << y_superior << endl;
        cout << "Y (inferior): " << setprecision(7) << y_inferior << endl;
        cout << "Error (superior): " << setprecision(7) << error_superior << endl;
        cout << "Error (inferior): " << setprecision(7) << error_inferior << endl;

        // Mostrar la y interpolada
        cout << "\ny = " << setprecision(7) << y_int << endl;

        // Mostrar el mejor error de interpolacion
        if (abs(error_superior) < abs(error_inferior)){
            cout << "\nMejor error de interpolacion: " << error_superior << endl;
        } else {
            cout << "\nMejor error de interpolacion: " << error_inferior << endl;
        }
    } 
}