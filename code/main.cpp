#include <iostream>
#include <vector>
#include "lagrange.h"
#include "newton.h"
#include "util.h"
#include "spline3.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;

using interpolacion::newton;
using interpolacion::lagrange;
using util::imprimir_tabla;
using interpolacion::spline3;

/** @brief Caso 1 Newton - Diapositivas */
void caso_1_newton();

/** @brief Caso 1 Lagrange - Diapositivas */
void caso_1_lagrange();

/** @brief Caso 2 Lagrange - Grado n */
void caso_2_lagrange();

/** @brief Caso 1 Trazadores Cubicos */
void caso_1_spline3();

int main(){
    // Menu para comparar los metodos de interpolacion de Newton y Lagrange con el caso 1

    int opcion;

    do{
        cout << "\n ------ Metodos de Interpolacion ------\n" << endl;
        cout << "1. Caso 1 Newton" << endl;
        cout << "2. Caso 1 Lagrange" << endl;
        cout << "3. Caso 2 Lagrange" << endl;
        cout << "4. Caso 1 Trazadores Cubicos" << endl;
        cout << "0. Salir" << endl;
        cout << "\nIngrese una opcion: ";
        cin >> opcion;

        switch(opcion){
            case 1:
                caso_1_newton();
                break;
            case 2:
                caso_1_lagrange();
                break;
            case 3:
                caso_2_lagrange();
                break;
            case 4:
                caso_1_spline3();
                break;
            case 0:
                cout << "Saliendo..." << endl;
                break;
            default:
                cout << "Opcion invalida" << endl;
                break;
        }
    } while(opcion != 0);

    return 0;
}

void caso_1_newton(){

    vector <double> x = {
        100.0f, 200.0f, 300.0f, 400.0f, 500.0f
    };

    vector <double> y = {
        -160.0f, -35.0f, -4.2f, 9.0f, 16.9f
    };

    // Instancia de Newton
    newton n(x, y);

    // Imprimir el polinomio
    cout << "\nPolinomio interpolante: " << n.polinomio() << endl;

    // Valor a interpolar
    double x_int;

    cout << "\nInterpolacion por diferencias divididas de Newton" << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, "Temperatura(K)", "B (cm3/mol)");
    // Solicitar el valor a interpolar
    do{
        cout << "Ingrese el valor a interpolar: ";
        cin >> x_int;
    } while(x_int < x[0] || x_int > x[x.size() - 1]);

    // Interpolar el valor ingresado por el usuario
    double y_int = n.interpolar(x_int);

    cout << "\ny = " << setprecision(7) << y_int << endl;
    
}

void caso_1_lagrange(){

    vector <double> x = {
        100.0f, 200.0f, 300.0f, 400.0f, 500.0f
    };

    vector <double> y = {
        -160.0f, -35.0f, -4.2f, 9.0f, 16.9f
    };

    // Instancia de Lagrange
    lagrange l(x, y);

    // Imprimir el polinomio
    cout << "\nPolinomio interpolante: " << l.polinomio() << endl;

    // Valor a interpolar
    double x_int;

    cout << "\nInterpolacion mediante el metodo de Lagrange" << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, "Temperatura(K)", "B (cm3/mol)");
    // Solicitar el valor a interpolar
    do{
        cout << "Ingrese el valor a interpolar: ";
        cin >> x_int;
    } while(x_int < x[0] || x_int > x[x.size() - 1]);

    // Interpolar el valor ingresado por el usuario
    double y_int = l.interpolar(x_int);

    cout << "\ny = " << setprecision(7) << y_int << endl;
    
}

void caso_2_lagrange(){

    vector <double> x = {
        100.0f, 200.0f, 300.0f, 400.0f, 500.0f
    };

    vector <double> y = {
        -160.0f, -35.0f, -4.2f, 9.0f, 16.9f
    };

    // Instancia de Lagrange
    lagrange l(x, y);

    // Imprimir el polinomio
    cout << "\nPolinomio interpolante: " << l.polinomio() << endl;

    // Valor a interpolar
    double x_int;
    int grado;

    cout << "\nInterpolacion mediante el metodo de Lagrange" << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, "Temperatura(K)", "B (cm3/mol)");

    // Solicitar el valor a interpolar
    do{
        cout << "Ingrese el valor a interpolar: ";
        cin >> x_int;
    } while(x_int < x[0] || x_int > x[x.size() - 1]);

    // Solicitar el grado de interpolación
    do{
        cout << "Ingrese el grado de interpolacion: ";
        cin >> grado;
    } while(grado < 1 || grado > x.size());

    // Interpolar el valor ingresado por el usuario con el grado especificado
    double y_int = l.interpolar(x_int, grado);

    cout << "\ny = " << setprecision(7) << y_int << endl;

    // Imprimir el error

    double error_int = abs(l.calcular_error_interpolacion(x_int, grado));

    cout << "\nError de interpolacion: " << error_int << endl;
    
}

void caso_1_spline3(){
    
    vector <double> x = {
        3.0f, 4.5f, 7.0f, 9.0f
    };

    vector <double> y = {
        2.5f, 1.0f, 2.5f, 0.5f
    };

    // Instancia de Trazador Cubico
    spline3 s3(x, y);

    // Imprimir el polinomio
    // cout << "Polinomio interpolante: " << l.polinomio() << endl;

    // Valor a interpolar
    double x_int;

    cout << "\nInterpolacion mediante el metodo de Trazadores Cubicos" << endl;

    // Imprimir la tabla
    imprimir_tabla(x, y, "  x  ", "  f(x)  ");
    // Solicitar el valor a interpolar
    do{
        cout << "Ingrese el valor a interpolar: ";
        cin >> x_int;
    } while(x_int < x[0] || x_int > x[x.size() - 1]);

    // Interpolar el valor ingresado por el usuario
    double y_int = s3.interpolar(x_int);

    cout << "\ny = " << setprecision(7) << y_int << endl;
    
}



