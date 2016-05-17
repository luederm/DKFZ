# author: Matthew Lueder
# description: This program shows basic usage of the Rcpp, inline, and, rbenchmark packages

#install.packages('Rcpp')
require('Rcpp')

#install.packages('rbenchmark')
require('rbenchmark')

#install.packages('inline')
require('inline')

# ***** Using inline/cxxfunction *****
# ** Example 1: Multiple each element of one vector with each of another **
funSrc <- '
// Implicitly calls "as" function
Rcpp::NumericVector xa(a);
Rcpp::NumericVector xb(b);
int n_xa = xa.size(), n_xb = xb.size();

Rcpp::NumericVector xab(n_xa + n_xb - 1);
for ( int i = 0; i < n_xa; ++i )
{
  for ( int j = 0; j < n_xb; ++j )
  {
     xab[i + j] += xa[i] * xb[j];
  }
  return xab;
}
'

inlineFun <- cxxfunction(sig = signature(a="numeric", b="numeric"),
                         funSrc,
                         plugin = "Rcpp")
inlineFun(1:22, 4:45)


# ** Example 2: Print out hello world N times **
funSrc2 <- '
#include <iostream>

int xb = Rcpp::as<int>(b);

for (int a = 0; a < xb; ++a)
  std::cout << "hello world" << std::endl;
'

iFunc2 <- cxxfunction(sig = signature(b="numeric"), funSrc2, plugin = "Rcpp")

iFunc2(10)
# Giving it more than one value will cause an error
iFunc2(1:3)
# Giving a decimal will cause digits after decimal to be truncated off
iFunc2(2.88)


# ** Example 3: Creating R datatypes in c++ **
funSrc3 <- '
using namespace Rcpp;

/* 
  The following objects inherit the RObject class. The RObject class is a c++ representation of an R object
  and is essentially a pointer around a SEXP object (SEXP is only data member). The constructor protect the object 
  from garbage collection, and the destructor removes protection. It has member functions for querying properties, 
  management of attributes, and handling slots. 
*/
CharacterVector a = CharacterVector::create("foo", "bar");
NumericVector b = NumericVector::create(0.1, 1.2);
LogicalVector c = LogicalVector::create(true, false, true, true);
IntegerVector d = IntegerVector::create(1, 2, 3);

List z = List::create(a, b, c, d);

// wrap function implicitly called in return statement
return z;
'
iFunc3 <- cxxfunction(sig = signature(), funSrc3, plugin = "Rcpp")
iFunc3()


# ** Example 4: Clone **
funSrc4 = '
#include <iostream>

// Clone() allows you to make copies of RObjects
Rcpp::NumericVector a = Rcpp::NumericVector::create(0, 0, 0);
Rcpp::NumericVector b = Rcpp::clone(a);

b[1] = 9.8; // This will not change vector a

List z = List::create(a, b);
return z;
'
iFunc4 <- cxxfunction(sig = signature(), funSrc4, plugin = "Rcpp")
iFunc4()


# ** Example 5: Matrices **
funSrc5 = '
Rcpp::NumericMatrix mat(inputMat);
std::transform(mat.begin(), mat.end(), mat.begin(), ::sqrt);

return mat;
'
input = matrix(c(2, 4, 3, 1, 5, 7), nrow=3, ncol=2)
iFunc5 <- cxxfunction(sig = signature(inputMat = "numeric"), funSrc5, plugin = "Rcpp" )
input
iFunc5(input)
# Clone was not used so the R object passed in was directly modified
input


# ** Example 6: Logical vectors **
funSrc6 = '
Rcpp::LogicalVector v(6);
/*
  R vectors support missing data, so this functionality exists in Rcpp`s logical vector.
  This means there are tree possible values rather than 2. NaN, Inf, and NA all collapse into NA.
*/
v[0] = false;
v[1] = true;
v[3] = R_NaN;
v[4] = R_PosInf;
v[5] = NA_REAL;

return v;
'
iFunc6 <- cxxfunction(sig = signature(), funSrc6, plugin = "Rcpp")
iFunc6()


# ** Example 7: The 'Named' helper class **
funSrc7 = '
Rcpp::NumericVector vec = Rcpp::NumericVector::create(
                            Rcpp::Named("mean") = 1.23,
                            Rcpp::Named("dim") = 42,
                            Rcpp::Named("cnt") = 12
                          );
return vec;
'
iFunc7 <- cxxfunction(sig = signature(), funSrc7, plugin = "Rcpp")
iFunc7()


# ** Example 8: Shortening syntax **
funSrc8 = '
// namespace Rcpp is implicitly defined when using plugin = "Rcpp" of cxxfunction.
// Rcpp::Named() can be shortened to _[]
NumericVector vec = NumericVector::create(
                      _["mean"] = 1.23,
                      _["dim"] = 42,
                      _["cnt"] = 12
                    );
return vec;
'
iFunc8 <- cxxfunction(sig = signature(), funSrc8, plugin = "Rcpp")
iFunc8()


# ** Example 9: Lists (GenericVectors) **
funSrc9 = '
/*
  We have seen how to make lists using `create` in previous examples, so here we will 
  first reserve space for elements of the list then assign them with operator[].
  Doing it this way automatically names reserved elements as an integer, starting at 1
  and counting up.
*/
List reservedList(3);
reservedList[0] = CharacterVector::create("foo", "bar");
reservedList[1] = IntegerVector::create(1, 2, 3);
reservedList[2] = LogicalVector::create(true, false, true, true);

// We can also use the functions push_back and push_front
reservedList.push_back("Last_Element");
reservedList.push_front("First_Element");

return reservedList;
'
iFunc9 <- cxxfunction(sig=signature(), funSrc9, plugin = "Rcpp")
iFunc9()


# ** Example 10: Data Frames **
funSrc10 = '
IntegerVector v = IntegerVector::create(7,8,9);
std::vector<std::string> s(3);
s[0] = "x";
s[1] = "y";
s[2] = "z";

return DataFrame::create(_["a"]=v, _["b"]=s);
'
iFunc10 <- cxxfunction(signature(), funSrc10, plugin="Rcpp")
iFunc10()


# ** Example 11: Functions **
funSrc11 = '
// In this simple example, a R function is passed in as an argument and used by C++ code
Function sort(x);
return sort( y, Named("decreasing", true));
'
iFunc11 <- cxxfunction(signature(x="function", y="ANY"), funSrc11, plugin="Rcpp")
iFunc11(sort, sample(1:5, 10, TRUE))
iFunc11(sort, sample(LETTERS[1:5], 10, TRUE))
# We can use a different function too, "decreasing" argument will be ignored
iFunc11(mean, sample(1:5, 10, TRUE))


# ** Example 12: Direct access of R functions **
funSrc12 = '
RNGScope scp;
// Rcpp functions are constructed with name of the R function you would like to access
Rcpp::Function rt("rt");
return rt(5, 3);
'
iFunc12 <- cxxfunction(signature(), funSrc12, plugin = "Rcpp")
iFunc12()


# ** Example 13: Enviroment Class **
funcSrc13 = '
Environment stats("package:stats");
Function rnorm = stats["rnorm"];
return rnorm(10, Named("sd", 100.0));
'
iFunc13 <- cxxfunction(signature(), funcSrc13, plugin = "Rcpp")
iFunc13()


# ** Example 14: Global Environment **
funcSrc14 = '
Environment global = Environment::global_env();

// Access variable in global env
float var = global["aVar"];

std::map<std::string,std::string> map;
map["foo"] = "oof";
map["bar"] = "rab";

// Add variable to global env
global["y"] = map;

return wrap(var);
'
iFunc14 <- cxxfunction(signature(), funcSrc14, plugin = "Rcpp")
aVar = 12345
iFunc14()
y


# ** Example 15: Inline - verbose **
funSrc15 = '
#include <iostream>

int xb = Rcpp::as<int>(b);

for (int a = 0; a < xb; ++a)
std::cout << "hello world" << std::endl;
'
# Set verbose = T
iFunc15 <- cxxfunction(sig = signature(b="numeric"), funSrc15, plugin = "Rcpp", verbose = T)
# After running this command you will see:
#   1) Environment variables set by cxxfunction
#   2) Location of linked files
#   3) The source code of the program
#   4) The compilation argument


