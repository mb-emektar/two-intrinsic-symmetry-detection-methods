#pragma once

#include "Display.h"
#include <memory> 
#include <iostream>
#include <windows.h>
#include <string>
#pragma warning(disable: 4996)

void myProgram(char* charptr)
{
    std::string strFilePath;

    std::cout << "Please Select the Symmetry Detection Method: " << std::endl;
    int methodSelection = 0;
    while (true) {
        std::cout << "1: Global Intrinsic Symmetries of Shapes " << std::endl;
        std::cout << "2: Möbius Transformations For Global Intrinsic Symmetry Analysis" << std::endl;
        std::cin >> methodSelection;

        if (methodSelection == 1 || methodSelection == 2) {
            break;
        }

        std::cout << "Wrong input, please select again. " << std::endl;
    }

    std::cout << "Please select the input: " << std::endl;
    int selection = 0;
    while (true) {
        std::cout << "1: man0 " << std::endl;
        std::cout << "2: horse0 " << std::endl;
        std::cout << "3: 251 " << std::endl;
        std::cout << "4: 398 " << std::endl;
        std::cout << "5: 200 " << std::endl;
        std::cout << "6: 185 " << std::endl;
        std::cout << "7: man " << std::endl;
        std::cout << "8: armadillo " << std::endl;
        std::cout << "9: homer " << std::endl;
        std::cout << "10: suzanne " << std::endl;
        std::cout << "11: dragon " << std::endl;
        std::cout << "12: centaur " << std::endl;
        std::cout << "13: horse " << std::endl;
        std::cout << "14: cheburashka " << std::endl;
        std::cout << "15: beast " << std::endl;
        std::cout << "16: spot " << std::endl;
        std::cout << "17: 0001.null.0" << std::endl;
        std::cout << "18: 0001.isometry.1 " << std::endl;
        std::cout << "19: 0001.isometry.2 " << std::endl;
        std::cout << "20: 0001.isometry.3 " << std::endl;
        std::cout << "21: 0001.isometry.4 " << std::endl;
        std::cout << "22: 0001.isometry.5 " << std::endl;
        std::cout << "23: 0001.isometry.6 " << std::endl;
        std::cout << "24: 0001.isometry.7 " << std::endl;
        std::cout << "25: 0001.isometry.8 " << std::endl;
        std::cout << "26: 0001.isometry.9 " << std::endl;
        std::cout << "27: 0001.isometry.10 " << std::endl;
        std::cout << "28: 0001.isometry.11 " << std::endl;
        std::cout << "29: 0001.isometry.12 " << std::endl;
        std::cout << "30: 0001.isometry.13 " << std::endl;
        std::cout << "31: 0001.isometry.14 " << std::endl;
        std::cout << "32: 0001.isometry.15 " << std::endl;
        std::cin >> selection;

        if (selection == 1) {
            strFilePath = "inputs/man0.off";
            break;
        }
        else if (selection == 2) {
            strFilePath = "inputs/horse0.off";
            break;
        }
        else if (selection == 3) {
            strFilePath = "inputs/251.off";
            break;
        }
        else if (selection == 4) {
            strFilePath = "inputs/398.off";
            break;
        }
        else if (selection == 5) {
            strFilePath = "inputs/200.off";
            break;
        }
        else if (selection == 6) {
            strFilePath = "inputs/185.off";
            break;
        }
        else if (selection == 7) {
            strFilePath = "inputs/man.off";
            break;
        }
        else if (selection == 8) {
            strFilePath = "inputs/armadillo.off";
            break;
        }
        else if (selection == 9) {
            strFilePath = "inputs/homer.off";
            break;
        }
        else if (selection == 10) {
            strFilePath = "inputs/suzanne.off";
            break;
        }
        else if (selection == 11) {
            strFilePath = "inputs/dragon.off";
            break;
        }
        else if (selection == 12) {
            strFilePath = "inputs/centaur.off";
            break;
        }
        else if (selection == 13) {
            strFilePath = "inputs/horse.off";
            break;
        }
        else if (selection == 14) {
            strFilePath = "inputs/cheburashka.off";
            break;
        }
        else if (selection == 15) {
            strFilePath = "inputs/beast.off";
            break;
        }
        else if (selection == 16) {
            strFilePath = "inputs/spot.off";
            break;
        }
        else if (selection == 17) {
            strFilePath = "inputs/0001.null.0.off";
            break;
        }
        else if (selection == 18) {
            strFilePath = "inputs/0001.isometry.1.off";
            break;
        }
        else if (selection == 19) {
            strFilePath = "inputs/0001.isometry.2.off";
            break;
        }
        else if (selection == 20) {
            strFilePath = "inputs/0001.isometry.3.off";
            break;
        }
        else if (selection == 21) {
            strFilePath = "inputs/0001.isometry.4.off";
            break;
        }
        else if (selection == 22) {
            strFilePath = "inputs/0001.isometry.5.off";
            break;
        }
        else if (selection == 23) {
            strFilePath = "inputs/0001.isometry.6.off";
            break;
        }
        else if (selection == 24) {
            strFilePath = "inputs/0001.isometry.7.off";
            break;
        }
        else if (selection == 25) {
            strFilePath = "inputs/0001.isometry.8.off";
            break;
        }
        else if (selection == 26) {
            strFilePath = "inputs/0001.isometry.9.off";
            break;
        }
        else if (selection == 27) {
            strFilePath = "inputs/0001.isometry.10.off";
            break;
        }
        else if (selection == 28) {
            strFilePath = "inputs/0001.isometry.11.off";
            break;
        }
        else if (selection == 29) {
            strFilePath = "inputs/0001.isometry.12.off";
            break;
        }
        else if (selection == 30) {
            strFilePath = "inputs/0001.isometry.13.off";
            break;
        }
        else if (selection == 31) {
            strFilePath = "inputs/0001.isometry.14.off";
            break;
        }
        else if (selection == 32) {
            strFilePath = "inputs/0001.isometry.15.off";
            break;
        }

        std::cout << "Wrong input, please select again. " << std::endl;
    }

    char* filePath = new char[strFilePath.length() + 1];
    std::strcpy(filePath, strFilePath.c_str());

    auto dsp = std::make_unique<Display>(charptr, filePath);

    dsp->createMembers(methodSelection);
    //dsp->calculateIsocurves();

    if (methodSelection == 1)
    {
        dsp->calculateIntrinsicSymmetry();
    }
    else if (methodSelection == 2) {
        dsp->calculateSampling();
        dsp->calculateMobiusSymmetry();
    }

    dsp->displayMesh();
}

int main(int, char** argv)
{
    std::cout << "Program instance started." << std::endl;

    myProgram(argv[0]);
    system("cls");

    std::cout << "Program instance exiting and restarting..." << std::endl;

    // Relaunch the application
    STARTUPINFO si;
    PROCESS_INFORMATION pi;

    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    // Get the path to the current executable
    wchar_t  path[MAX_PATH];
    GetModuleFileName(NULL, path, MAX_PATH);

    if (CreateProcess(path, NULL, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) {
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);
    }
	return 0;
}