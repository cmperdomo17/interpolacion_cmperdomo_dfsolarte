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

using regresion::solucion_lineal;
using regresion::solucion_potencia;
using regresion::solucion_exponencial;
using regresion::solucion_cuadratica;

using regresion::lineal_simple;
using regresion::potencia;
using regresion::exponencial;
using regresion::cuadratica;


/** @brief Metodo Newton - Grado n */
void interpolacion_newton(vector<double> x, vector<double> y);

/** @brief Metodo Trazadores Cubicos */
void interpolacion_spline3(vector<double> x, vector<double> y);

/**
* @brief Metodo 1 Regresion Lineal Simple
* @param x Valores de la variable independiente
* @param y Valores de la variable dependiente
* @param title Titulo del caso
* @param x_label Etiqueta de la variable independiente para la tabla
* @param y_label Etiqueta de la variable dependiente para la tabla
*/
void regresion_lineal_simple(vector<double> x, 
							 vector<double> y, 
							 string title, 
							 string x_label="", 
							 string y_label="");

/**
* @brief Metodo 2 Regresion Linealizada mediante la funcion potencia
* @param x Valores de la variable independiente
* @param y Valores de la variable dependiente
* @param title Titulo del caso
* @param x_label Etiqueta de la variable independiente para la tabla
* @param y_label Etiqueta de la variable dependiente para la tabla
*/
void regresion_potencia(vector<double> x, 
							 vector<double> y, 
							 string title, 
							 string x_label="", 
							 string y_label="");

/**
* @brief Metodo 3 Regresion Linealizada mediante la funcion exponencial
* @param x Valores de la variable independiente
* @param y Valores de la variable dependiente
* @param title Titulo del caso
* @param x_label Etiqueta de la variable independiente para la tabla
* @param y_label Etiqueta de la variable dependiente para la tabla
*/
void regresion_exponencial(vector<double> x, 
						   vector<double> y, 
						   string title, 
						   string x_label="", 
						   string y_label="");

/**
* @brief Metodo 4 Regresion cuadratica
* @param x Valores de la variable independiente
* @param y Valores de la variable dependiente
* @param title Titulo del caso
* @param x_label Etiqueta de la variable independiente para la tabla
* @param y_label Etiqueta de la variable dependiente para la tabla
*/
void regresion_cuadratica(vector<double> x, 
						  vector<double> y, 
						  string title, 
						  string x_label="", 
						  string y_label="");

int main(){
    // Menu para comparar los metodos de interpolacion de Newton y Lagrange con el caso 1
    /*
        Vectores temporales -- Borrar
    */
    vector <double> x = {0.4f, 0.8f, 1.3f, 1.8f, 2.0f, 2.2f, 2.6f};
    vector <double> y = {1.452360f, 1.995632f, 2.719678f, 3.273019f, 3.359425f, 3.316678f, 2.669452f};

    int opcion;

    do{
        cout << "\n ------ Metodos de Interpolacion ------\n" << endl;
        cout << "1. Metodo Newton" << endl;
        cout << "2. Metodo Trazadores Cubicos" << endl;
        cout << " ------ Metodos de Regresion ------\n" << endl;
        cout << "3. Metodo Lineal Simple" << endl;
        cout << "4. Metodo Potencial" << endl;
        cout << "5. Metodo Exponencial" << endl;
        cout << "6. Metodo Cuadratico" << endl;
        cout << "0. Salir" << endl;
        cout << "\nIngrese una opcion: ";
        cin >> opcion;

        switch(opcion){
            case 1:
                interpolacion_newton(x,y);
                break;
            case 2:
                interpolacion_spline3(x,y);
                break;
            case 3:
                regresion_lineal_simple(x,y,"3");
                break;
            case 4:
                regresion_potencia(x,y,"4");
                break;
            case 5:
                regresion_exponencial(x,y,"5");
                break;
            case 6:
                regresion_cuadratica(x,y,"6");
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

void interpolacion_newton(vector<double> x, vector<double> y){
    
        // Instancia de Newton
        newton n(x, y);
    
        // Imprimir el polinomio
        cout << "\nPolinomio interpolante: " << n.polinomio() << endl;
    
        // Valor a interpolar
        double x_int;
        size_t grado;
    
        cout << "\nInterpolacion por diferencias divididas de Newton" << endl;
    
        // Imprimir la tabla
        imprimir_tabla(x, y, "Temperatura(K)", "B (cm3/mol)");
    
        // Solicitar el valor a interpolar
        do{
            cout << "Ingrese el valor a interpolar: ";
            cin >> x_int;
        } while(x_int < x[0] || x_int > x[x.size() - 1]);
    
        // Solicitar el grado de interpolación
        do{
            cout << "Ingrese el grado de interpolacion (0 para utilizar todos los datos): ";
            cin >> grado;
        } while(grado < 0 || grado > x.size());
        
        // Interpolar el valor ingresado por el usuario con el grado especificado
        double y_int;

        if(grado==0){
            y_int= n.interpolar(x_int);
        }else{
            y_int = n.interpolar(x_int, grado);
        }
    
        n.interpolar(x_int, grado);
    
        cout << "\ny = " << setprecision(7) << y_int << endl;
    
        // Imprimir el error
    
        double error_int = abs(n.calcular_error_interpolacion(x_int, grado));
    
        cout << "\nError de interpolacion: " << error_int << endl;
}

void interpolacion_spline3(vector<double> x, vector<double> y){
    
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

void regresion_lineal_simple(vector<double> x, 
							 vector<double> y, 
							 string title, 
							 string x_label, 
							 string y_label){
	cout << title << endl; 
	
	//Imprimir tabla
	imprimir_tabla(x,y,x_label,y_label);
	
	//Crear una instancia de regresion simple
	lineal_simple ls(x,y);
	solucion_lineal sol = ls.calcular();
	
	//Imprimir la solucion
	sol.imprimir();
}
	
void regresion_potencia(vector<double> x, 
							 vector<double> y, 
							 string title, 
							 string x_label, 
							 string y_label){
	cout << title << endl; 
	
	//Imprimir tabla
	imprimir_tabla(x,y,x_label,y_label);
	
	//Crear una instancia de regresion simple
	potencia reg_potencia(x,y);
	
	//Calcular la solucion
	solucion_potencia sol = reg_potencia.calcular();
	
	//Imprimir la solucion
	sol.imprimir();
}
	
void regresion_exponencial(vector<double> x, 
						vector<double> y, 
						string title, 
						string x_label, 
						string y_label){
	cout << title << endl; 
	
	//Imprimir tabla
	imprimir_tabla(x,y,x_label,y_label);
	
	//Crear una instancia de regresion simple
	exponencial reg_exponencial(x,y);
	
	//Calcular la solucion
	solucion_exponencial sol = reg_exponencial.calcular();
	
	//Imprimir la solucion
	sol.imprimir();
}

void regresion_cuadratica(vector<double> x, 
						   vector<double> y, 
						   string title, 
						   string x_label, 
						   string y_label){
	cout << title << endl; 
	
	//Imprimir tabla
	imprimir_tabla(x,y,x_label,y_label);
	
	//Crear una instancia de regresion cuadratica
	cuadratica reg_cuadratica(x,y);
	
	//Calcular la solucion
	solucion_cuadratica sol = reg_cuadratica.calcular();
	
	//Imprimir la solucion
	sol.imprimir();
}