#include "vtkLAADFilter.h"

#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkStaticPointLocator.h>
#include <vtkInformationDoubleKey.h>

#include <cmath>

#include <vtkSmartPointer.h>
#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkStandardNewMacro(vtkLAADFilter);

vtkLAADFilter::vtkLAADFilter() {}
vtkLAADFilter::~vtkLAADFilter() {}

int vtkLAADFilter::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}

int vtkLAADFilter::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}

void vtkLAADFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

int vtkLAADFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *outInfo2 = outputVector->GetInformationObject(1);

  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  this->computeLAAD(input, output);

  return 1;
}

void vtkLAADFilter::computeLAAD(vtkUnstructuredGrid *input, vtkUnstructuredGrid *output)
{
    vtkIdType numPoints = input->GetNumberOfPoints();

    // Pass on data to output
    output->CopyStructure(input);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetCellData()->PassData(input->GetCellData());

    // Create new array
    VTK_CREATE(vtkDoubleArray, newData);
    newData->SetName("LAAD");
    newData->SetNumberOfComponents(1);
    newData->SetNumberOfTuples(numPoints);

    // Get pointers to the velocity component arrays
    vtkDoubleArray *globalU = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetAbstractArray("u"));
    vtkDoubleArray *globalV = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetAbstractArray("v"));
    vtkDoubleArray *globalW = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetAbstractArray("w"));

    VTK_CREATE(vtkStaticPointLocator, locator);
    locator->SetDataSet(input);
    locator->BuildLocator();

    // According to the VTK documentation, the GetPoint(vtkIdType, double x[3]) method
    // is only thread safe when first called from a single thread. GetPoint is invoked
    // here solely for that reason, and serves no other purpose.
    double tmp[3];
    input->GetPoint(0, tmp);

    // Main loop over all points in the dataset
    #pragma omp parallel for
    for(vtkIdType pointId = 0; pointId < numPoints; pointId++)
    {
        double thisCoords[3];
        input->GetPoint(pointId, thisCoords);

        VTK_CREATE(vtkIdList, neighbourhood);
        locator->FindPointsWithinRadius(this->Radius, thisCoords, neighbourhood);
        vtkIdType neighbourhoodSize = neighbourhood->GetNumberOfIds();

        double uRef = 0;
        double vRef = 0;
        double wRef = 0;

        // Compute average vector in neighbourhood
        for (vtkIdType i = 0; i < neighbourhoodSize; i++)
        {
            vtkIdType id = neighbourhood->GetId(i);
            uRef += globalU->GetValue(id);
            vRef += globalV->GetValue(id);
            wRef += globalW->GetValue(id);
        }

        uRef /= neighbourhoodSize;
        vRef /= neighbourhoodSize;
        wRef /= neighbourhoodSize;

        double normRef = std::sqrt(uRef * uRef + vRef * vRef + wRef * wRef);

        double LAADValueAtIndex = 0;
        double uOther, vOther, wOther, normOther, cosAngle;

        for (vtkIdType i = 0; i < neighbourhoodSize; i++)
        {
            vtkIdType id = neighbourhood->GetId(i);
            uOther = globalU->GetValue(id);
            vOther = globalV->GetValue(id);
            wOther = globalW->GetValue(id);

            normOther = std::sqrt(uOther * uOther + vOther * vOther + wOther * wOther);

            if (normRef * normOther == 0) continue;

            cosAngle = (uRef * uOther + vRef * vOther + wRef * wOther) / (normRef * normOther);
            if (cosAngle < -1) cosAngle = -1;
            if (cosAngle > 1) cosAngle = 1;

            LAADValueAtIndex += std::acos(cosAngle);
        }

        LAADValueAtIndex /= neighbourhoodSize * M_PI;

        newData->SetValue(pointId, LAADValueAtIndex);
    }

    output->GetPointData()->AddArray(newData);

}
