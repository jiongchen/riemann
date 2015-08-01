#include "ipopt_solver.h"

#include <coin/IpNLP.hpp>

namespace riemann {

class nlp : public Ipopt::NLP
{

};

}
