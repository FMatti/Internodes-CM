template <typename T>
class PenaltyFunctionLinear {
public:
  T operator()(T& var) { return var; }
  T derivative(T& var) { return 1.0; }
};
template <typename T>
class PenaltyFunctionQuadratic {
public:
  T operator()(T& var) { return var*var + var; }
  T derivative(T& var) { return 2.0*var + 1.0; }
};