#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <cmath>
#include <iostream>
#include <vtkType.h>

int main() {
    //  Parametre kruhu
    const int nPoints = 100;        // počet bodov na kružnici
    const double cx = 0.5;          // stred X
    const double cy = 0.5;          // stred Y
    const double r  = 0.3;          // polomer

    //  Vytvorenie bodov na kružnici
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < nPoints; ++i) {
        double theta = 2.0 * M_PI * i / nPoints;
        double x = cx + r * std::cos(theta);
        double y = cy + r * std::sin(theta);
        double z = 0.0;  // 2D krivka leží v rovine XY
        points->InsertNextPoint(x, y, z);
    }

    //  Polyline (spojenie bodov do uzavretej krivky)
    vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
    polyLine->GetPointIds()->SetNumberOfIds(nPoints + 1); // +1, aby sa uzavrel okruh
    for (int i = 0; i < nPoints; ++i)
        polyLine->GetPointIds()->SetId(i, i);
    polyLine->GetPointIds()->SetId(nPoints, 0); // späť na prvý bod

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(polyLine);

    //  PolyData objekt obsahujúci krivku
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(cells);

    //  Zápis do .vtp súboru
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("circle.vtp");
    writer->SetInputData(polyData);
    writer->Write();
    std::cout << " Soubor 'circle.vtp' bol vytvorený.\n";

    //****************************
    // Vytvorenie kartezskej siete
    vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();
    int extent = 10; 
    grid->SetExtent(0, extent, 0, extent, 0, 0);
    grid->SetOrigin(0.0, 0.0, 0.0);
    double dx = 0.1, dy = 0.1;
    grid->SetSpacing(dx, dy, 0.1);

    // Alokácia pamäte pre skalárne dáta
    grid->AllocateScalars(VTK_DOUBLE, 1); 

    // Celkový počet bodov
    const int nx = extent + 1;
    const int ny = extent + 1;

    // Nastavenie všetkých hodnôt na 0
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double* ptr = static_cast<double*>(grid->GetScalarPointer(i, j, 0));
            *ptr = 0.0;
        }
    }

    // Zápis siete do .vti súboru
    vtkSmartPointer<vtkXMLImageDataWriter> gridWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    gridWriter->SetFileName("cartesian_grid.vti");
    gridWriter->SetInputData(grid);
    gridWriter->Write();
    std::cout << "Soubor 'cartesian_grid.vti' bol vytvorený.\n";

    return 0;
}