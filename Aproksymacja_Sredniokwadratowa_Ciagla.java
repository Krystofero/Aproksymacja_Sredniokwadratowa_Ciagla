package Aproksymacja_Sredniokwadratowa_Ciagla;

import java.util.Scanner;

public class Aproksymacja_Sredniokwadratowa_Ciagla {

    public static double funkcja(double x) {
        // return (double) Math.sqrt(x);
        return (double) (Math.pow(x,3) + 2);
    }
    public static double[] macierz_gauss(double[][] A, double[] b) {
        int N  = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double   t    = b[p]; b[p] = b[max]; b[max] = t;

            if (Math.abs(A[p][p]) <= 1e-10) {
                throw new RuntimeException("Macierz jest pojedyncza lub prawie pojedyncza");
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }
//Kwadratura Gaussa-Legendre'a -----------------------------------------------
    static protected double[] w30 = {0.1028526528935588, //wagi kwadratury
            0.1028526528935588,
            0.1017623897484055,
            0.1017623897484055,
            0.0995934205867953,
            0.0995934205867953,
            0.0963687371746443,
            0.0963687371746443,
            0.0921225222377861,
            0.0921225222377861,
            0.0868997872010830,
            0.0868997872010830,
            0.0807558952294202,
            0.0807558952294202,
            0.0737559747377052,
            0.0737559747377052,
            0.0659742298821805,
            0.0659742298821805,
            0.0574931562176191,
            0.0574931562176191,
            0.0484026728305941,
            0.0484026728305941,
            0.0387991925696271,
            0.0387991925696271,
            0.0287847078833234,
            0.0287847078833234,
            0.0184664683110910,
            0.0184664683110910,
            0.0079681924961666,
            0.0079681924961666};

    static protected double[] t30 = {-0.0514718425553177, //węzły kwadratury
            0.0514718425553177,
            -0.1538699136085835,
            0.1538699136085835,
            -0.2546369261678899,
            0.2546369261678899,
            -0.3527047255308781,
            0.3527047255308781,
            -0.4470337695380892,
            0.4470337695380892,
            -0.5366241481420199,
            0.5366241481420199,
            -0.6205261829892429,
            0.6205261829892429,
            -0.6978504947933158,
            0.6978504947933158,
            -0.7677774321048262,
            0.7677774321048262,
            -0.8295657623827684,
            0.8295657623827684,
            -0.8825605357920527,
            0.8825605357920527,
            -0.9262000474292743,
            0.9262000474292743,
            -0.9600218649683075,
            0.9600218649683075,
            -0.9836681232797472,
            0.9836681232797472,
            -0.9968934840746495,
            0.9968934840746495};

    public static double calka1(double a, double b, int stopien1, int stopien2){

        double h, h2;
        h = (b - a) / 2;
        h2 = (b + a) / 2;
        double wynik = 0;
        int stopien;
        stopien = stopien1 + stopien2;//mamy dwie funkcje bazowe(ten sam x do jakiejś potęgi), dlatego możemy dodać stopnie potęg do siebie
                                    //x.Math.pow(z)*x.Math.pow(z)=x.Math.pow(z+z)
        for (int i = 0; i < 30; i++) {
            wynik += Math.pow(((h) * t30[i] + h2),stopien) * w30[i];
        }
        return wynik*h;
    }

    public static double calka2(double a, double b,int stopien){

        double h, h2;
        h = (b - a) / 2;
        h2 = (b + a) / 2;
        double wynik = 0;

        for (int i = 0; i < 30; i++) {
            wynik += Math.pow(((h) * t30[i] + h2),stopien)* funkcja(((h) * t30[i] + h2)) * w30[i];
        }

        return wynik*h;
    }

    public static void main(String[] args) {
        double a, b; //przedział w jakim aproksymjemy funkcję
        double x; // x, w którym poszukujemy rozwiązania
        int st; //stopień wielomianu aproksymującego

        Scanner sc = new Scanner(System.in);
        System.out.println("Podaj dolny przedzial calkowania (a): ");
        a = Integer.parseInt(sc.nextLine());
        System.out.println("Podaj górny przedzial calkowania (b): ");
        b = Integer.parseInt(sc.nextLine());
        System.out.println("Podaj liczbe st (stopien wielomianu): ");
        st = Integer.parseInt(sc.nextLine());
//        st = 2;
        System.out.println("Podaj liczbę x, aby policzyć dla niej wartość funkcji aproksymującej");
        x = Double.parseDouble(sc.nextLine());
        sc.close();

        double[][] gauss=new double[st+1][st+2];
        double[] gauss2=new double[st+1];
        for (int i = 0; i < st+1; i++) {   //wyliczamy całki i otrzymujemy układ równań z niewiadomymi a
            for (int j = 0; j < st+1; j++) {
                gauss[i][j]=calka1(a,b,i,j);
                System.out.print(gauss[i][j]+" ");
            }
            gauss2[i]=calka2(a,b,i);
            System.out.println(" | " + gauss2[i]);
        }

        double[] wspolczynnik = new double[st+1]; //współczynniki a
        wspolczynnik=macierz_gauss(gauss,gauss2); //Eliminacją Gaussa wyliczam współczynniki

        for (int i = 0; i < st+1; i++) {
            System.out.println("wspolczynniki (a) "+wspolczynnik[i]);
        }

        System.out.println("Wynik funkcji(bez aproksymacji): "+ funkcja(x));

        double wynik=0;

        for (int i=0; i < st+1;i++) //współczynniki mnożę razy funkcje bazowe i wszystko sumuję
        {
            wynik+=wspolczynnik[i]*Math.pow(x,i);
        }
        System.out.println("Wynik po aproksymacji: "+wynik);
    }
}