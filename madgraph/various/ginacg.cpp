#include <ginac/ginac.h>
using namespace GiNaC;
#include <iostream>
using namespace std;
#include <sstream>
using namespace std;
#include <string>

extern "C" {

float GinacG(float y[], float xx,int w)
{
    int i;
    lst l = {};
    ex x= numeric(xx);
    for(i=0;i<w;i++)
        l.append(numeric(y[i]));
    ostringstream s;
    s<< G(l,x).evalf();

    return stof(s.str());
}
}