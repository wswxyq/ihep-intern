/*
Writing a new p.d.f class using RooClassFactory
The utility class RooClassFactory simplifies the task of writing a custom RooFit p.d.f class to writing
the actual p.d.f expression.
The class factory has several modes of operation. The simplest mode of operation is for a function
expression that is simple enough to be expressed in a single line of code. For those cases, a
completely functional custom p.d.f. class can be written as follows: 
*/

RooClassFactory::makePdf("RooMyPdf","x,alpha",0,"sqrt(abs(x*alpha))+0.1");

/*
This operation writes out two files, RooMyPdf.cxx and RooMyPdf.h, that can be compiled and linked
with ACliC for immediate use

root>.L RooMyPdf.cxx+ 

If your function expression is not simple enough to be expressed in a single line of code, you can
simply omit the expression when you request RooClassFactory to create the class 
*/
//no space betwwen items!!
RooClassFactory::makePdf("RooMyPdf","x,alpha") ; right!
RooClassFactory::makePdf("RooMyPdf","x, alpha") ; wrong!

/*
This creates a fully functional class with a dummy implementation of RooAbsPdf::evaluate(). To
make it a functional class, edit the file RooMyPdf.cxx and insert a function expression as return value
in the evaluate() method of your class, using as many lines of code as you need. 
*/

//remember to change return value!!!

Double_t RooMyPdf::evaluate() const
 {
 // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
 return 1 ;
 } 

/*
This creates a fully functional class with a dummy implementation of RooAbsPdf::evaluate(). To
make it a functional class, edit the file RooMyPdf.cxx and insert a function expression as return value
in the evaluate() method of your class, using as many lines of code as you need. 

Since RooAbsPdf have no fixed interpretation of variables being observables or parameters, there is
no need, or point, in explicitly normalizing the expression in evaluate() with respect to a specific
choice of observables: the return value of evaluate() is always divided a posteriori by a
normalization integral before it is return through RooAbsPdf::getVal().

By default this normalization step is done using a numeric integrator, but if you know how to integrate
your class over one (or more) choices of observables, you can advertise this capability in the p.d.f.
and your analytical integral will be used instead of numeric integration whenever it is appropriate. You
can invoke RooClassFactory::makePdf() with different options that will make skeleton code for
the analytical integral interface. Details can be found in the RooClassFactory HTML class
documentation. 
*/