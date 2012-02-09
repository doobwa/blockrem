# I have a class that stores data and various statistics about my data.  I want to be able to work with it from R and from C++ functions, but I do not want to have the entire data structure created as an R object.  I thought one might be able to return a pointer from the object to R, and pass that pointer from R to a C++ function that needs to interact with the object.  This seems to work until garbage collection.  

# I have included an example below illustrating my current approach.  I altered the World example to return a pointer to the current object.  I have a C++ function that makes changes to the object.

# I saw this post http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2011-December/003214.html  but I wanted to use modules and I haven't been able to adapt that solution to my situation.



library(inline)
library(Rcpp)
fx <- cxxfunction(,"",includes=
  '
#include <Rcpp.h>

class World {
public:
  World() : msg("hello") {}
  void set(std::string msg) { 
    this->msg = msg; 
  }
  std::string greet() { 
    return msg; 
  }
  SEXP ptr() {
   // return wrap(XPtr<World>(this, true));
    return XPtr<World>(this, false);
  }

private:
  std::string msg;
};

int fn(SEXP ptr_) {
  World *s = XPtr<World>(ptr_);
  s->set("c++ function has been here");
  return 1;
}

RCPP_MODULE(example){
using namespace Rcpp ;
function("fn", &fn);
class_<World>( "World")
  .constructor()
  .method( "greet", &World::greet , 
  "get the message")
  .method( "set", &World::set , 
  "set the message")
  .method( "ptr", &World::ptr , 
  "get a pointer");
}
', plugin="Rcpp")

example <- Module("example",getDynLib(fx))

s <- new(example$World)

# Interact with World object from R
s$greet()
s$set("hello from R")
s$greet()

# Grab pointer to this World object
s$ptr()  

# Call a c++ function that uses World object
example$fn(s$ptr())
s$greet()  # c++ function has altered s, as desired

# Causes segfault
gc()


fx <- cxxfunction(,"",includes=
  '
using namespace Rcpp;
template<class T>
  XPtr<T> unwrap_robject(const SEXP& sexp){
    RObject ro(sexp);
    if(ro.isObject()){
      Language call("as.environment",sexp);
      SEXP ev = call.eval();
      Language call1("get",".pointer",-1,ev);
      SEXP ev1 = call1.eval();
      XPtr<T > xp(ev1);
      return xp;
    }else{
      XPtr<T > xp(sexp);
      return xp;
    }
  }

template<class T>
  SEXP wrap_in_reference_class(const T& obj,std::string class_name){
    XPtr< T > xp(new T(obj));
    Language call( "new", Symbol( class_name ),xp);
    return call.eval();
  }

class MyObject{
  public:
    MyObject(SEXP sexp){
      XPtr<MyObject> xp = unwrap_robject<MyObject>(sexp);
      /* copy members from xp here */
    }
  
  operator SEXP() const{
    return wrap_in_reference_class(*this,"MyObject");
  }
};

  RCPP_MODULE(MyModule){
    
    class_<MyObject>("MyObject")
    .constructor<SEXP>();
    
  }
', plugin="Rcpp")

example <- Module("MyModule",getDynLib(fx))
m <- new(example$MyObject)