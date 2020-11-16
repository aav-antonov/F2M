#ifndef Ax_b_h
#define AX_b_h
//////////////////////

#include <unordered_map>



class Ax_b {

    private:


    std::vector<std::vector<int> > Ai;
    std::vector<std::vector<int> > Av;
    std::vector<int> b  ;
    std::vector<double> d, Aj;
    std::vector<double> x_grad , b_grad;

    std::unordered_map<int,int> x_stat;

    

    double ax(int j , std::vector<double> &X);
    double product(std::vector<double> &X1, std::vector<double> &X2);
    double product(std::vector<int> &X1, std::vector<double> &X2);
    void x_grad_update();
    void b_grad_update();
    void x_update();
    void d_j(int j);
    void d_update();

    void print(std::vector<double> &X1);
    void vector2null(std::vector<double> &X1);

    void parse_Ax_b(std::vector<int> & A);

   public:
    double d_max_error , d_max_error_STOP = 1e-8;
    int go_max = 0; 
    int j_count = 0, i_count = 0;
    std::vector<double> x;
    void solve();
    void print();
   

    //Ax_b(std::vector<std::vector<int>> &Ai , std::vector<std::vector<int>> &Av, std::vector<int> &b){Ai = Ai; Av = Av;b = b;};
    Ax_b(std::string file);
    Ax_b(std::vector<int> & A){parse_Ax_b(A);};
    ~Ax_b(){}//std::cout << "destructur Ax_b\n";

};








    

#endif

