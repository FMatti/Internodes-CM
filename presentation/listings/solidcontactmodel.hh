class SolidContactModel : public Model {
public:
  virtual ~SolidContactModel() = 0;

  // Polymorphic accessors
  virtual Model & getContactMechanicsModel() = 0;
  virtual AbstractContactDetector & getContactDetector() = 0;
  virtual SolidMechanicsModel & getSolidMechanicsModel() = 0;

  // Helper methods, for example:
  template <typename FunctorType>
  inline void applyBC(const FunctorType & func) {
    getSolidMechanicsModel().applyBC(func);
  }
  // ...
};