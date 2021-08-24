#include "Params.h"
#include <unordered_map>
#include <string>

Params::Params() {
    intParamMap.clear();
    doubleParamMap.clear();
}

// Params::Params(const Params &s) {
//     this->intParamMap = s.intParamMap;
//     this->doubleParamMap = s.doubleParamMap;
// }

void Params::addParam(string &paramName, double val) {
    doubleParamMap[paramName] = val;
}

void Params::addParam(string &paramName, int val) {
    // doubleParamMap[paramName] = val;
    intParamMap[paramName] = val;
}

void Params::getParam(string &paramName, int &paramVal) {
    paramVal = intParamMap[paramName];
}

void Params::getParam(string &paramName, double &paramVal) {
    paramVal = doubleParamMap[paramName];
}

void Params::addParam(const char *paramName, double val) {
    string temp(paramName);
    this->addParam(temp, val);

}

void Params::addParam(const char *paramName, int val) {
    string temp(paramName);
    this->addParam(temp, val);

}

void Params::getParam(const char *paramName, double &paramVal) {
    string temp(paramName);
    this->getParam(temp, paramVal);

}

void Params::getParam(const char *paramName, int &paramVal) {
    string temp(paramName);
    this->getParam(temp, paramVal);
}

