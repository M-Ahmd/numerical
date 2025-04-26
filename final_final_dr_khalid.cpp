#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

enum FunctionType {
    Power,
    Log,
    Exp,
    Sqrt,
    Sin,
    Cos,
    Tan,
    Sec,
    Csc,
    Cot,
    Asin,
    Acos,
    Atan,
    Sqr,
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Sech,
    Csch,
    Asinh,
    Acosh,
    Atanh,
    Acoth,
    Asech,
    Acsch
};

double acoth(double x) {
    return 0.5 * log((x + 1) / (x - 1));
}

double asech(double x) {
    return log(1 / x + sqrt(1 / (x * x) - 1));
}

double acsch(double x) {
    return log(1 / x + sqrt(1 / (x * x) + 1));
}

struct Term;

struct Factor {
    FunctionType type;
    int power;
    vector<Term> innerExpr;
    
    bool isComposite()const
    {
        return !innerExpr.empty();
    }
};

struct Term {
    double coeff;
    vector<Factor> factors;

    double eval(double x) const {
        double result = coeff;
        for (int i = 0; i < factors.size(); ++i) {
            double factorValue = 1.0;
            
            if (factors[i].isComposite())
            {
                
                double innerValue = 0.0;
                vector<Term> composite_terms = factors[i].innerExpr;
                for (int j = 0; j < composite_terms.size(); j++) {
                    innerValue += composite_terms[j].eval(x);
                }
                
                
                switch (factors[i].type) {
                    case Power:   factorValue = pow(innerValue, factors[i].power); break;
                    case Log:     factorValue = log(innerValue); break;
                    case Exp:     factorValue = exp(innerValue); break;
                    case Sqrt:    factorValue = sqrt(innerValue); break;
                    case Sin:     factorValue = sin(innerValue); break;
                    case Cos:     factorValue = cos(innerValue); break;
                    case Tan:     factorValue = tan(innerValue); break;
                    case Sec:     factorValue = 1.0 / cos(innerValue); break;
                    case Csc:     factorValue = 1.0 / sin(innerValue); break;
                    case Cot:     factorValue = 1.0 / tan(innerValue); break;
                    case Asin:    factorValue = asin(innerValue); break;
                    case Acos:    factorValue = acos(innerValue); break;
                    case Atan:    factorValue = atan(innerValue); break;
                    case Sqr:     factorValue = pow(innerValue, 2); break;
                    case Sinh:    factorValue = sinh(innerValue); break;
                    case Cosh:    factorValue = cosh(innerValue); break;
                    case Tanh:    factorValue = tanh(innerValue); break;
                    case Coth:    factorValue = 1.0 / tanh(innerValue); break;
                    case Sech:    factorValue = 1.0 / cosh(innerValue); break;
                    case Csch:    factorValue = 1.0 / sinh(innerValue); break;
                    case Asinh:   factorValue = asinh(innerValue); break;
                    case Acosh:   factorValue = acosh(innerValue); break;
                    case Atanh:   factorValue = atanh(innerValue); break;
                    case Acoth:   factorValue = acoth(innerValue); break;
                    case Asech:   factorValue = asech(innerValue); break;
                    case Acsch:   factorValue = acsch(innerValue); break;
                    default: break;
                }
            } else {
                
                switch (factors[i].type) {
                    case Power:   factorValue = pow(x, factors[i].power); break;
                    case Log:     factorValue = log(x); break;
                    case Exp:     factorValue = exp(x); break;
                    case Sqrt:    factorValue = sqrt(x); break;
                    case Sin:     factorValue = sin(x); break;
                    case Cos:     factorValue = cos(x); break;
                    case Tan:     factorValue = tan(x); break;
                    case Sec:     factorValue = 1.0 / cos(x); break;
                    case Csc:     factorValue = 1.0 / sin(x); break;
                    case Cot:     factorValue = 1.0 / tan(x); break;
                    case Asin:    factorValue = asin(x); break;
                    case Acos:    factorValue = acos(x); break;
                    case Atan:    factorValue = atan(x); break;
                    case Sqr:     factorValue = pow(x, 2); break;
                    case Sinh:    factorValue = sinh(x); break;
                    case Cosh:    factorValue = cosh(x); break;
                    case Tanh:    factorValue = tanh(x); break;
                    case Coth:    factorValue = 1.0 / tanh(x); break;
                    case Sech:    factorValue = 1.0 / cosh(x); break;
                    case Csch:    factorValue = 1.0 / sinh(x); break;
                    case Asinh:   factorValue = asinh(x); break;
                    case Acosh:   factorValue = acosh(x); break;
                    case Atanh:   factorValue = atanh(x); break;
                    case Acoth:   factorValue = acoth(x); break;
                    case Asech:   factorValue = asech(x); break;
                    case Acsch:   factorValue = acsch(x); break;
                    default: break;
                }
            }
            if(factors[i].type != Power)
                result *= pow(factorValue, factors[i].power);
            else
                result *= factorValue;
        }
        return result;
    }
};

vector<Term> createExpression() {
    vector<Term> expr;
    int T;
    cout << "Enter number of terms: ";
    cin >> T;

    while (T-- > 0) {
        Term t;
        cout << "Enter coefficient: ";
        cin >> t.coeff;

        int F;
        cout << "Enter number of factors: ";
        cin >> F;

        while (F-- > 0) {
            Factor f;
            cout << "Choose type:\n";
            cout << "1=x^n, 2=ln, 3=exp, 4=sqrt, 5=sin, 6=cos, 7=tan,\n";
            cout << "8=sec, 9=csc, 10=cot, 11=arcsin, 12=arccos, 13=arctan,\n";
            cout << "14=sqr(x), 15=sinh(x), 16=cosh(x), 17=tanh(x), 18=csch(x),\n";
            cout << "19=sech(x), 20=coth(x), 21=arcsinh(x), 22=arccosh(x),\n";
            cout << "23=arctanh(x), 24=arcoth(x), 25=arsech(x), 26=arcsch(x)\n";
            cout << "Enter type: ";
            int typeInput;
            cin >> typeInput;

            if (typeInput >= 1 && typeInput <= 26) {
                f.type = static_cast<FunctionType>(typeInput - 1);
            } else {
                cout << "Invalid input! Please enter a valid function type.\n";
                continue;
            }

            cout << "Enter power: ";
            cin >> f.power;

            char composite;
            cout << "Is this a composite function? (y/n): ";
            cin >> composite;

            if (composite == 'y' || composite == 'Y') {
                cout << "Creating inner expression for this function...\n";
                f.innerExpr = createExpression();
            }
            t.factors.push_back(f);
        }
        expr.push_back(t);
    }
    return expr;
}

double evalExpr(double x, const vector<Term>& expr) {
    double sum = 0;
    for (int i = 0; i < expr.size(); ++i) {
        sum += expr[i].eval(x);
    }
    return sum;
}

double bisectionWithTable(double a, double b, double tol, const vector<Term>& expr, int& nIter) {
    double fa = evalExpr(a, expr);
    double fb = evalExpr(b, expr);

    if (fa * fb >= 0) {
        cout << "Error: f(a) * f(b) must be negative for the bisection method.\n";
        return -1;
    }

    double mid;

    cout << setw(5) << "Iter" << setw(15) << "a" << setw(15) << "b" << setw(15) << "p_n" << setw(15) << "f(p_n)" << endl;

    for(int i = 1; i <= nIter; i++) {
        mid = (a + b) / 2;
        double fm = evalExpr(mid, expr);

        cout << setw(5) << i << setw(15) << fixed << setprecision(9) << a
             << setw(15) << b << setw(15) << mid << setw(15) << fm << endl;

        if (fabs(fm) < tol) {
            nIter = i;
            break;
        }

        if (fa * fm < 0) {
            b = mid;
            fb = fm;
        } else {
            a = mid;
            fa = fm;
        }
    }

    return mid;
}

int main() {
    cout << "Building your function...\n";
    vector<Term> expr = createExpression();

    
    double testX;
    cout << "Enter a value to test your function: ";
    cin >> testX;
    cout << "f(" << testX << ") = " << evalExpr(testX, expr) << endl;

    
    double a, b, tol;
    cout << "Enter interval a and b (a < b): ";
    cin >> a >> b;

    if (a >= b) {
        cout << "Error: Ensure that a < b for the interval.\n";
        return 1;
    }

    cout << "Enter tolerance: ";
    cin >> tol;

    int nIter = ceil((log(b - a) - log(tol)) / log(2));
    double root = bisectionWithTable(a, b, tol, expr, nIter);

    
    cout << fixed << setprecision(9);
    cout << "Approximate root = " << root << endl;
    cout << "Number of iterations = " << nIter << endl;
    

    return 0;
}