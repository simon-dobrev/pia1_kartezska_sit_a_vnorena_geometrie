#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>

#include <vtkSelectEnclosedPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolygon.h>

#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vtkType.h>

int main() {
    
    //  Parametre kruhu
    const int nPoints = 100;        // počet bodov na kružnici
    const double cx = 0.5;          // stred X
    const double cy = 0.5;          // stred Y
    const double r  = 0.3;          // polomer

    const double TOLERANCE = 1e-3;  // tolerancia pre testovanie vnútra kruhu

    //  Vytvorenie bodov na kružnici
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < nPoints; ++i) {
        double theta = 2.0 * M_PI * i / nPoints;
        double x = cx + r * std::cos(theta);
        double y = cy + r * std::sin(theta);
        double z = 0.0;  // 2D krivka leží v rovine XY
        points->InsertNextPoint(x, y, z);
    }
/*
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
*/

    //Pouziti Polygon misto PolyLine

    	vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(nPoints);
	for (int i = 0; i < nPoints; ++i)
    		polygon->GetPointIds()->SetId(i, i);

	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polygon);

	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(points);
	polyData->SetPolys(cells);

    //  Zápis do .vtp súboru
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("circle.vtp");
    writer->SetInputData(polyData);
    writer->Write();
    std::cout << " Soubor 'circle.vtp' bol vytvorený.\n";

    //****************************
    // Vytvorenie kartezskej siete
    vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();
    int extent = 100; 
    grid->SetExtent(0, extent, 0, extent, 0, 0);
    grid->SetOrigin(0.0, 0.0, 0.0);
    double dx = 1.0 / extent, dy = 1.0 / extent;
    grid->SetSpacing(dx, dy, 0.1);

    // Celkový počet bodov
    const int nx = extent + 1;
    const int ny = extent + 1;

    // Nastavenie všetkých hodnôt na 0
    /*for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double* ptr = static_cast<double*>(grid->GetScalarPointer(i, j, 0));
            
	    int ijk[3] = {i, j, 0}; // místo {i, j, 0}
	    
	    vtkIdType pid = grid->ComputePointId({ijk});
            int inside = selector->IsInside(pid);
	    *ptr = inside ? 1.0 : 0.0; // 1 = uvnitř kruhu, 0 = venku
        }
    }*/

	// Pole pro cell data
vtkSmartPointer<vtkDoubleArray> cellValues = vtkSmartPointer<vtkDoubleArray>::New();
cellValues->SetName("InsideCircle");
cellValues->SetNumberOfTuples(grid->GetNumberOfCells());
	
	// Super-sampling: míra pokrytí buňky kruhem
    
    int samples = 10; // jemnost mřížky uvnitř buňky
    std::ofstream fout("gridCoverage.txt");

    for (int j = 0; j < extent; ++j) {
        for (int i = 0; i < extent; ++i) {
            double x0 = i * dx;
            double x1 = (i + 1) * dx;
            double y0 = j * dy;
            double y1 = (j + 1) * dy;

            int insideCount = 0;
            int totalCount = samples * samples;

            for (int sj = 0; sj < samples; ++sj) {
                for (int si = 0; si < samples; ++si) {
                    double xs = x0 + (si + 0.5) * (dx / samples);
                    double ys = y0 + (sj + 0.5) * (dy / samples);

                    double dx_c = xs - cx;
                    double dy_c = ys - cy;
                    /*if (dx_c*dx_c + dy_c*dy_c <= r*r) {
                        insideCount++;
                    }*/
                    // Vypočítajte druhú mocninu vzdialenosti od stredu
                    double distSq = dx_c * dx_c + dy_c * dy_c;

                    // Testujeme, či je rozdiel vzdialenosti a polomeru menší ako tolerancia.
                    if (std::abs(distSq - r*r) <= TOLERANCE) { 
                        insideCount++;
                    }
                }
            }

        double coverage = static_cast<double>(insideCount) / totalCount;
       	vtkIdType cellId = j * extent + i;
      	cellValues->SetValue(cellId, coverage);	

	fout << static_cast<int>(coverage * 100) << "%\t";
    }
	fout << "\n";
}
	fout.close();
    	grid->GetCellData()->AddArray(cellValues);

    // Zápis siete do .vti súboru
    vtkSmartPointer<vtkXMLImageDataWriter> gridWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    gridWriter->SetFileName("cartesian_grid.vti");
    gridWriter->SetInputData(grid);
    gridWriter->Write();
    std::cout << "Soubor 'cartesian_grid.vti' bol vytvorený.\n";

    return 0;
}
