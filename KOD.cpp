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
    
    // Parametre kruhu
    const int nPoints = 100;
    const double cx = 0.5;
    const double cy = 0.5;
    const double r  = 0.3;

    // -------- Kruh: body + polygon --------
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < nPoints; ++i) {
        double theta = 2.0 * M_PI * i / nPoints;
        double x = cx + r * std::cos(theta);
        double y = cy + r * std::sin(theta);
        double z = 0.0;
        points->InsertNextPoint(x, y, z);
    }

    // Vytvoření polygonu
    vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
    polygon->GetPointIds()->SetNumberOfIds(nPoints);
    for (int i = 0; i < nPoints; ++i)
        polygon->GetPointIds()->SetId(i, i);

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(polygon);

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(cells);

    // Zápis kruhu do .vtp súboru
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("circle.vtp");
    writer->SetDataModeToBinary();
    writer->SetInputData(polyData);
    writer->Write();

    // -------- Karterzska mriežka --------
    vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();
    int extent = 100; // Počet buněk v jednom směru
    
    // Extent v VTK je definován přes body (n buněk = n+1 bodů)
    grid->SetExtent(0, extent, 0, extent, 0, 0);
    grid->SetOrigin(0.0, 0.0, 0.0);

    double dx = 1.0 / extent;
    double dy = 1.0 / extent;
    grid->SetSpacing(dx, dy, 0.1);

    vtkIdType nCells = grid->GetNumberOfCells();

    // Cell data: pokrytie
    vtkSmartPointer<vtkDoubleArray> coverage =
        vtkSmartPointer<vtkDoubleArray>::New();
    coverage->SetName("InsideCircleCoverage");
    coverage->SetNumberOfTuples(nCells);

    // Cell data: normála
    vtkSmartPointer<vtkDoubleArray> normals =
        vtkSmartPointer<vtkDoubleArray>::New();
    normals->SetName("NormalVector");
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(nCells);

    // Inicializácia normál
    double zero[3] = {0.0, 0.0, 0.0};
    for (vtkIdType i = 0; i < nCells; ++i)
        normals->SetTuple(i, zero);

    // -------- Super-sampling --------
    int samples = 10;
    int total = samples * samples;
    std::ofstream fout("gridCoverage.txt");

    for (int j = 0; j < extent; ++j) {
        for (int i = 0; i < extent; ++i) {

            double x0 = i * dx;
            double y0 = j * dy;

            int inside = 0;

            // Super-sampling cyklus
            for (int sj = 0; sj < samples; ++sj) {
                for (int si = 0; si < samples; ++si) {
                    double xs = x0 + (si + 0.5) * (dx / samples);
                    double ys = y0 + (sj + 0.5) * (dy / samples);

                    double dx_c = xs - cx;
                    double dy_c = ys - cy;
		     
                    if (dx_c * dx_c + dy_c * dy_c <= r * r)
                        inside++;
                }
            }

            double cov = static_cast<double>(inside) / total;
            vtkIdType cellId = j * extent + i;
            coverage->SetValue(cellId, cov);
        }
    }

    // Přidání pole pokrytí do gridu
    grid->GetCellData()->AddArray(coverage);

    // -------- Získání normál před gradient --------
    // Normála se počítá jako záporný gradient pole pokrytí: n = -grad(coverage)
    // Iterujeme o 1 méně na krajích, abychom mohli sahat na sousedy (i-1, i+1)

    for (int j = 1; j < extent - 1; ++j) {
        for (int i = 1; i < extent - 1; ++i) {
            
            vtkIdType cellId = j * extent + i;
            double cov = coverage->GetValue(cellId);

            // Normály počítáme pouze na rozhraní dle pokrytí
            if (cov > 0.0 && cov < 1.0) {
                
                // Indexy sousedních buněk (Left, Right, Down, Up)
                vtkIdType id_L = j * extent + (i - 1);
                vtkIdType id_R = j * extent + (i + 1);
                vtkIdType id_D = (j - 1) * extent + i;
                vtkIdType id_U = (j + 1) * extent + i;

                double val_L = coverage->GetValue(id_L);
                double val_R = coverage->GetValue(id_R);
                double val_D = coverage->GetValue(id_D);
                double val_U = coverage->GetValue(id_U);

                // Gradient míří dovnitř, proto použijeme záporné znaménko 
                double grad_x = (val_R - val_L) / (2.0 * dx);
                double grad_y = (val_U - val_D) / (2.0 * dy);

                double nx = -grad_x;
                double ny = -grad_y;

                // Normalizace vektoru
                double mag = std::sqrt(nx * nx + ny * ny);
                if (mag > 1e-9) {
                    nx /= mag;
                    ny /= mag;
                } else {
                    // Pokud je gradient nulový, normálu nenastavujeme
                    nx = 0.0; ny = 0.0;
                }

                double n[3] = {nx, ny, 0.0};
                normals->SetTuple(cellId, n);
            }
        }
    }

    grid->GetCellData()->AddArray(normals);

    // Zápis výsledkov do .vti súboru
    vtkSmartPointer<vtkXMLImageDataWriter> gridWriter =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    gridWriter->SetFileName("cartesian_grid.vti");
    gridWriter->SetDataModeToBinary();
    gridWriter->SetInputData(grid);
    gridWriter->Write();

    std::cout << "Hotovo. Vytvoren soubor 'cartesian_grid.vti'" << std::endl;

    return 0;
}
