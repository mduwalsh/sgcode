#include <cstdio>
#include <iostream>

using namespace std;

int main()
{
  cout << "Hello World!"<<endl;
  char myntcs[] = "some text";
string mystring = myntcs;  // convert c-string to string
cout << mystring <<endl;          // printed as a library string
cout << mystring.c_str() <<endl;  // printed as a c-string 
sprintf(myntcs,"wer");
cout << myntcs << endl;
cout << mystring.c_str() <<endl;
int x;
int y = 10;
const int * p = &y;
x = *p;          // ok: reading p
y = 2;
cout<<x<<endl<<y<<endl<<*p<<endl;     
char str[300];
sprintf(str, "%d%.2f%s", 3, 4.566, "random string");
string ss = str;
cout<<str<<endl;
}
