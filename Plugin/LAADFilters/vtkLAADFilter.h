#ifndef __vtkLAADFilter_h
#define __vtkLAADFilter_h

#include <vtkUnstructuredGridAlgorithm.h>
#include "vtkSystemIncludes.h"
#include <set>

class vtkLAADFilter : public vtkUnstructuredGridAlgorithm
{
public:
    static vtkLAADFilter *New();
    // Sets up IsTypeOf, IsA, SafeDownCast
    vtkTypeMacro(vtkLAADFilter,vtkUnstructuredGridAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent) override;

    // Radius of the spherical neighbourhood
    double Radius = 0;

    // Macros to get/set the radius, and rerun computations if modified
    vtkSetMacro(Radius, double);
    vtkGetMacro(Radius, double);

protected:
  vtkLAADFilter();
  ~vtkLAADFilter();

  virtual int FillInputPortInformation(int port, vtkInformation *info) override;
  virtual int FillOutputPortInformation(int port, vtkInformation* info) override;

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override;

  void computeLAAD(vtkUnstructuredGrid *, vtkUnstructuredGrid *);
  void computeNeighbourhood(vtkUnstructuredGrid *, double, vtkIdType, std::set<vtkIdType> *);

private:
  vtkLAADFilter(const vtkLAADFilter&);
  void operator=(const vtkLAADFilter&);
};

#endif
