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
        double z = 0.0;  // 2D krivka leží v rovine XY
        points->InsertNextPoint(x, y, z);
    }

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
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("circle.vtp");
    writer->SetDataModeToBinary();
    writer->SetInputData(polyData);
    writer->Write();

    // -------- Karterzska mriežka --------
    vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();
    int extent = 100;
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

            for (int sj = 0; sj < samples; ++sj)
                for (int si = 0; si < samples; ++si) {
                    double xs = x0 + (si + 0.5) * (dx / samples);
                    double ys = y0 + (sj + 0.5) * (dy / samples);

                    double dx_c = xs - cx;
                    double dy_c = ys - cy;

                    if (dx_c * dx_c + dy_c * dy_c <= r * r)
                        inside++;
                }

            double cov = static_cast<double>(inside) / total;
            vtkIdType cellId = j * extent + i;
            coverage->SetValue(cellId, cov);

            double nx = 0.0, ny = 0.0;

            if (cov > 0.0 && cov < 1.0) {
                nx = (x0 + 0.5 * dx) - cx;
                ny = (y0 + 0.5 * dy) - cy;

                double mag = sqrt(nx * nx + ny * ny);
                if (mag > 1e-9) {
                    nx /= mag;
                    ny /= mag;
                }

                double n[3] = {nx, ny, 0.0};
                normals->SetTuple(cellId, n);
            }

            fout << static_cast<int>(cov * 100) << "%\t"
                 << nx << "\t" << ny << "\t";
        }
        fout << "\n";
    }
    fout.close();

    grid->GetCellData()->AddArray(coverage);
    grid->GetCellData()->AddArray(normals);

    // Zápis siete do .vti súboru
    vtkSmartPointer<vtkXMLImageDataWriter> gridWriter =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    gridWriter->SetFileName("cartesian_grid_norm.vti");
    gridWriter->SetDataModeToBinary();
    gridWriter->SetInputData(grid);
    gridWriter->Write();

    return 0;
}
