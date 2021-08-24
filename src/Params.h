#ifndef PARAMS_H
#define PARAMS_H

#include <unordered_map>
#include <string>

using namespace std;

class Params {
public:
    Params();
    void addParam(string &paramName, double val);
    void addParam(string &paramName, int val);

    void addParam(const char *paramName, double val);
    void addParam(const char *paramName, int val);

    void getParam(string &paramName, double &paramVal);
    void getParam(string &paramName, int &paramVal);

    void getParam(const char *paramName, double &paramVal);
    void getParam(const char *paramName, int &paramVal);
protected:
    unordered_map<string, double> doubleParamMap;
    unordered_map<string, int> intParamMap;
};

#endif