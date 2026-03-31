#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>

using namespace std;

ifstream fin("EVOLUTIE.IN");
ofstream fout("EVOLUTIE.out");

int n, a, b, coef_A, coef_B, coef_C, prob_recombinare, prob_mutatie, etape, precizie;

int nr_biti;
double pas_discretizare;
double probabilitatea[1005], capat_selectie[1005];
bool prima_generatie = true;

struct Individ{
    string cromozom;
    double x;
    double fitness;
};


Individ populatie[1005];

void citire(){
    cout<<"Dimensiunea populatiei: ";
    fin>>n;
    cout<<"Coeficientul 'a' al ecuatiei a * x^2 + b * x + c : ";
    fin>>coef_A;
    cout<<"Coeficientul 'b' al ecuatiei "<<coef_A<<" * x^2 + b * x + c : ";
    fin>>coef_B;
    cout<<"Coeficientul 'c' al ecuatiei "<<coef_A<<" x^2 + "<<coef_B<<" * x + c : ";
    fin>>coef_C;
    cout<<"Functia : "<<coef_A<<" * x^2 + "<<coef_B<<" * x + "<<coef_C<<'\n';
    cout<<"Precizie : ";
    fin>>precizie;
    cout<<"Domeniul capat stanga [ _ , ] : ";
    fin>>a;
    cout<<"Domeniul capat dreapta:["<<a<<" , _ ] : ";
    fin>>b;
    cout<<"Domeniul : [ "<<a<<" , "<<b<<" ]"<<'\n';
    cout<<"Probabilitatea de recombinare (? %) : ";
    fin>>prob_recombinare;
    cout<<"Probabilitatea de  mutatie (? %) : ";
    fin>>prob_mutatie;
    cout<<"Numar etape : ";
    fin>>etape;
}

string random_bits(){
    string s(nr_biti, '0');
    for (char &c : s){
        c = '0' + rand() % 2;
    }
    return s;
}

double decode(string cromozom){
    long long k = 0;
    for (char c : cromozom)
        k = k * 2 + (c - '0');
    return a + k * pas_discretizare;
}

double fitness(double x){
    return coef_A * x * x + coef_B * x + coef_C;
}

void afiseaza_populatie_initiala() {
    if (prima_generatie == true){
        fout<<"Populatia initiala"<<'\n';
        for (int i = 1; i <= n; i++) {
            fout<<setw(3)<<i<<" : "<<populatie[i].cromozom<<"  x = "<<setw(10)<<populatie[i].x<<"  f = "<<setw(10)<<populatie[i].fitness<<'\n';
        }
    }
}

void selectie(){
    double sum = 0;
    for (int i=1; i<=n; i++)
        sum += populatie[i].fitness;

    for (int i=1; i<=n; i++)
        probabilitatea[i] = populatie[i].fitness / sum;

    if (prima_generatie == true){
        fout<<'\n';
        fout<<"Probabilitati selectie"<<'\n';
        for (int i=1; i<=n; i++)
            fout<<setw(11)<<"cromozom "<<setw(4)<<i<<" probabilitate "<<setw(12)<<probabilitatea[i]<<'\n';
    }
}

void intervale_selectie(){
    for (int i=1; i<=n; i++){
        capat_selectie[i] = i==1 ? 0 : capat_selectie[i-1] + probabilitatea[i-1];
    }

    if (prima_generatie == true){
        fout<<'\n';
        fout<<"Intervale probabilitati selectie"<<'\n';
        for (int i=1; i<n; i++){
            fout<<"[ "<<capat_selectie[i]<<" , "<<capat_selectie[i+1]<<" ]"<<'\n';
        }
        fout<<"[ "<<capat_selectie[n]<<" , "<<"1.0"<<" ]"<<'\n';
    }
}

int cauta_interval(double u) {
    int stanga = 1, dreapta = n;
    while (stanga < dreapta) {
        int mid = (stanga + dreapta) / 2;
        if (u < capat_selectie[mid+1]) dreapta = mid;
        else stanga = mid + 1;
    }
    return stanga;
}

void selecteaza_populatie_noua() {
    Individ populatie_noua[1005];
    
    int elite=1;
    for (int i=2; i<=n; i++)
        if (populatie[i].fitness > populatie[elite].fitness) 
            elite=i;
    
    if (prima_generatie == true) 
        fout<<'\n';
    for (int i=1; i<=n; i++){
        double random = (double)rand() / RAND_MAX;
        int index_cromozom = cauta_interval(random);

        if (prima_generatie == true)
            fout<<"u = "<<setprecision(16)<<random<< "  ->  selectam cromozomul "<<index_cromozom<<'\n';
        populatie_noua[i] = populatie[index_cromozom];
    }
    populatie_noua[1] = populatie[elite];

    for (int i = 1; i <= n; i++)
        populatie[i] = populatie_noua[i];

    if (prima_generatie == true){
        fout<<'\n';
        fout<<"Dupa selectie:\n";
        for (int i = 1; i <= n; i++)
            fout<<setw(3)<<i<<" : "<<populatie[i].cromozom<<"  x = "<<setw(10)<<setprecision(6)<<populatie[i].x<<"  f = "<<setw(10)<<setprecision(6)<<populatie[i].fitness<<'\n';
    }
}

void crossover(){
    double prob_cross = (double)prob_recombinare/100;
    if (prima_generatie == true){
        fout<<'\n';
        fout<<"Probabilitatea de incrucisare "<<prob_cross<<'\n';
    }

    int participare[1005], k=0;
    if (prima_generatie == true)
        fout<<'\n';
    for (int i=1; i<=n; i++){
        double random = (double)rand() / RAND_MAX;

        if (prima_generatie == true){
            fout<<setw(3)<<i<<": " <<populatie[i].cromozom<< "  u = " <<setprecision(16)<<random;
            if (random < prob_cross) fout<<" < "<<prob_cross<<"  participa";
            fout<<'\n';
        }

        if (random < prob_cross) participare[++k] = i;

    }

    for (int j=1; j+1<=k; j+=2){
        int x = participare[j];
        int y = participare[j+1];

        int punct_rupere = 1 + rand() % (nr_biti - 1);

        if (prima_generatie == true){
            fout<<"\nRecombinare intre cromozomul "<<x<< " si cromozomul "<<y<<":\n";
            fout<<populatie[x].cromozom<<"  "<<populatie[y].cromozom<<"  punct = "<<punct_rupere<<'\n';
        }

        string newX = populatie[x].cromozom.substr(0, punct_rupere) + populatie[y].cromozom.substr(punct_rupere);
        string newY = populatie[y].cromozom.substr(0, punct_rupere) + populatie[x].cromozom.substr(punct_rupere);

        populatie[x].cromozom = newX;
        populatie[x].x = decode(newX);
        populatie[x].fitness = fitness(populatie[x].x);

        populatie[y].cromozom = newY;
        populatie[y].x = decode(newY);
        populatie[y].fitness = fitness(populatie[y].x);

        if (prima_generatie == true)
            fout<<"Rezultat:  "<<populatie[x].cromozom<<"  "<<populatie[y].cromozom<<'\n';
    }

    if (prima_generatie == true){
        fout<<'\n';
        fout<<"Dupa recombinare"<<'\n';
        for (int i=1;i<=n;i++){
            fout<<setw(3)<<i<<" : "<<populatie[i].cromozom<<"  x = "<<setw(10)<<setprecision(6)<<populatie[i].x<<"  f = "<<setw(10)<<setprecision(6)<<populatie[i].fitness<<'\n';
        }
    }
}

void mutatie(){
    double prob_mut = (double)prob_mutatie / 100;

    if (prima_generatie == true){
        fout<<'\n';
        fout<<"Probabilitate de mutatie pentru fiecare gena "<<prob_mut<<'\n';
    }

    Individ elita = populatie[1];

    if (prima_generatie == true)
        fout<<"Au fost modificati cromozomii :"<<'\n';
    for (int i=2; i<=n; i++){
        bool ok = false;
         for (char& bit : populatie[i].cromozom) {
            double random = (double)rand() / RAND_MAX;
            if (random < prob_mut){ 
                bit ^= 1; 
                ok = true;
            }
        }
        if (ok == true){
            populatie[i].x = decode(populatie[i].cromozom);
            populatie[i].fitness = fitness(populatie[i].x);
            if (prima_generatie == true)
                fout<<i<<'\n';
        }
    }
    populatie[1] = elita;

    
    if (prima_generatie == true){
        fout << "\nDupa mutatie:\n";
        for (int i = 1; i <= n; i++){
            fout<<setw(3)<<i<<" : "<<populatie[i].cromozom<<"  x = "<<setw(10)<<setprecision(6)<<populatie[i].x<<"  f = "<<setw(10)<<setprecision(6)<<populatie[i].fitness<<'\n';
        }
    }
}

int main(){
    citire();

    nr_biti = ceil(log2((b - a) * pow(10, precizie)));
    pas_discretizare = (b-a)/pow(2,nr_biti);

    for (int i=1; i<=n; i++){
        populatie[i].cromozom =  random_bits();
        populatie[i].x = decode(populatie[i].cromozom);
        populatie[i].fitness = fitness(populatie[i].x);
    }

    afiseaza_populatie_initiala();
    selectie();
    intervale_selectie();
    selecteaza_populatie_noua();
    crossover();
    mutatie();

    prima_generatie = false;

    fout<<'\n';
    fout<<"Evolutia maximului"<<'\n';
    for (int generatie=2; generatie<=etape; generatie++){
        afiseaza_populatie_initiala();
        selectie();
        intervale_selectie();
        selecteaza_populatie_noua();
        crossover();
        mutatie();

        double Max = populatie[1].fitness, medie = 0;
        for (int i=1; i<=n; i++) {
            Max = max(Max, populatie[i].fitness);
            medie += populatie[i].fitness;
        }
        medie /= n;
        fout<<"Generatia "<<setw(3)<<generatie<<" :  Max Fitness = "<<setprecision(15)<<Max<<"  Mean Fitness = "<<medie<<'\n';
    }
}
