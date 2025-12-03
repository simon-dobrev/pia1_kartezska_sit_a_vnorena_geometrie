# Nastavenie minimálnej verzie CMake
cmake_minimum_required(VERSION 3.12)

# Názov projektu
project(pia1_program LANGUAGES CXX)

# Nastavenie štandardu C++ na modernú verziu (odporúčané C++14 alebo 17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Zapnutie upozornení kompilátora
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# =================================================================
# 1. Nájdenie VTK
# 
# Vyžadujeme špecifické VTK moduly. Je kľúčové, aby ste tu uviedli 
# VŠETKY moduly, ktoré používate vo vašom C++ kóde.
# =================================================================
find_package(VTK REQUIRED COMPONENTS
    # Pre PolyData, Points, CellArray
    vtkCommonCore
    # Pre PolyDataNormals (Linker sa sťažoval, že táto chýba)
    vtkFiltersCore 
    # Pre Polygón a dáta buniek
    vtkCommonDataModel 
    # Pre vtkImageData
    vtkImagingCore
    # Pre vtkXMLPolyDataWriter a vtkXMLImageDataWriter
    vtkIOXML
)

# Ak sa VTK nenájde, vypíše chybu
if (NOT VTK_FOUND)
  message(FATAL_ERROR "Nebola nájdená knižnica VTK.")
endif()

# =================================================================
# 2. Definícia spustiteľného súboru
# =================================================================
# Nahraďte 'KOD.cpp' názvom vášho zdrojového súboru
add_executable(pia1_program KOD.cpp)

# =================================================================
# 3. Prepojenie knižníc (Linkovanie)
# =================================================================
target_link_libraries(pia1_program 
    PUBLIC 
        # Pripojenie všetkých nájdených VTK knižníc
        ${VTK_LIBRARIES}
)

# Pridanie hlavičkových súborov VTK
target_include_directories(pia1_program
    PRIVATE
        ${VTK_INCLUDE_DIRS}
)