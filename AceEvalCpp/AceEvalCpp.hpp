#pragma once

#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

using std::vector;
using std::string;
using std::set;
using std::map;
using std::istream;
using std::ifstream;
using std::exception;
using std::invalid_argument;
using std::runtime_error;


class Variable{
private:
  string fName;
  vector<string> fDomainNames;
  
public: //formerly protected
  Variable(string name, vector<string> domainNames){
    fName = name;
    fDomainNames.assign(domainNames.begin(), domainNames.end());
  }
  
public:
  string name() const {
    return fName;
  }
  
  vector<string> domainNames() const {
    return fDomainNames;
  }
  
  //TODO overload << instead
  string toString() {
    return name();
  }
  
  // functions to be used by std::map

  bool operator< (const Variable& other) const {
    return this->fName < other.name();
  }
  
  Variable(){
  }
  
  Variable(const Variable& other){
    fName = other.name();
    fDomainNames = other.domainNames();
  }
};

class Potential {
private:
    string fName;
    int fNumPositions;
    
public: //formerly protected
    Potential(string name, int numPositions) {
        fName = name;
        fNumPositions = numPositions;
    }
    
public:
    string name() const{
        return fName;
    }
    
    int numPositions() const{
        return fNumPositions;
    }
    
    //TODO overload << instead
    string toString() {
        return "T (" + name() + ")";
    }
    
    
    // functions to be used by std::map
    
    bool operator<(const Potential& other) const {
      return this->fName < other.name();
    }
    
    Potential(){
    }
    
    Potential(const Potential& other){
      fName = other.name();
      fNumPositions = other.numPositions();
    }
    
};

class OnlineEngine;

class Evidence {
  
private:
  OnlineEngine& fEngine;
  
public: //formerly protected:
  vector<double> fVarToCurrentNegWeight;
  vector<double> fVarToCurrentPosWeight;
  
private:
  double defaultWeight (int l);
  void setCurrentWeight (int l, double w);
  void setCurrentWeightToDefault (int l);
  void setCurrentWeightsToDefaults (const vector<int>&);
  void setCurrentWeights (double w, const vector<int>&);
  
public:
  void retractAll ();
  Evidence (OnlineEngine& engine);
  void varCommit (Variable v, int u);
  void valCommit (Variable v, int u, bool allow);
  void varRetract (Variable v);
  void varSet (Variable v, double w);
  void parmCommit (Potential t, int p, double w);
  void parmRetract (Potential t, int p);
  
};

class OnlineEngine {
public: // formerly protected:
    static const char CONSTANT = 0;
    static const char LITERAL = 1;
    static const char MULTIPLY = 2;
    static const char ADD = 3;
    vector<char> fNodeToType;
    vector<int> fNodeToLastEdge;
    vector<int> fNodeToLit;
    vector<int> fEdgeToTailNode;
    vector<int> fVarToNegLitNode;
    vector<int> fVarToPosLitNode;
    static const string READ_DELIMITER;
    static const string DELIMITER;
    set<Variable> fVariables;
    set<Potential> fPotentials;
    map<string, Variable> fNameToSrcVar;
    map<string, Potential> fNameToSrcPot;
    vector<double> fLogicVarToDefaultNegWeight;
    vector<double> fLogicVarToDefaultPosWeight;
    map<Variable, vector<int>* > fSrcVarToSrcValToIndicator;
    map<Potential, vector<int>* > fSrcPotToSrcPosToParameter;
    vector<double> fNodeToValue;
    vector<double> fNodeToDerivative;
    vector<bool> fNodeToOneZero;
    bool fUpwardPassCompleted;
    bool fTwoPassesCompleted;
    vector<double> fAcVarToMostRecentNegWeight;
    vector<double> fAcVarToMostRecentPosWeight;

    vector<double> clone(vector<double> a);
    void readArithmeticCircuit(istream& r);
    void readLiteralMap(istream& r);
    void upwardPass (Evidence ev);
    void twoPasses (Evidence ev);
    double computedValue (int n);
    int first (int n);
    int rootNode ();
    int numAcNodes ();

public:    
    OnlineEngine (string acFilename, string lmFilename);
    void initialize(istream& acReader, istream& lmReader);
    Variable varForName (string n);
    Potential potForName (string n);
    set<Variable> variables ();
    set<Potential> potentials ();
    void assertEvidence (Evidence e, bool secondPass);
    double probOfEvidence ();
    vector<double> varPartials (Variable v);
    map<Variable,vector<double> > varPartials (set<Variable> vs);
    vector<double> varMarginals (Variable v);
    map<Variable,vector<double> > varMarginals (set<Variable> vs);
    vector<double> varPosteriors (Variable v);
    map<Variable,vector<double> > varPosteriors (set<Variable> vs);
    vector<double> potPartials (Potential pot);
    map<Potential,vector<double> > potPartials (set<Potential> ps);
    vector<double> potMarginals (Potential p);
    map<Potential,vector<double> > potMarginals (set<Potential> ps);
    vector<double> potPosteriors (Potential p);
    map<Potential,vector<double> > potPosteriors (set<Potential> ps);
};

double Evidence::defaultWeight (int l) {
    return l < 0 ? fEngine.fLogicVarToDefaultNegWeight[-l]:
           fEngine.fLogicVarToDefaultPosWeight[l];
}

void Evidence::setCurrentWeight (int l, double w) {
    if (l < 0) {
        fVarToCurrentNegWeight[-l] = w;
    } else {
        fVarToCurrentPosWeight[l] = w;
    }
}

void Evidence::setCurrentWeightToDefault (int l) {
    setCurrentWeight (l, defaultWeight (l));
}

void Evidence::setCurrentWeightsToDefaults (const vector<int>& lits) {
    for (int l : lits) {
        setCurrentWeightToDefault (l);
    }
}

void Evidence::setCurrentWeights (double w, const vector<int>& lits) {
    for (int l : lits) {
        setCurrentWeight (l, w);
    }
}



void Evidence::retractAll () {
    fVarToCurrentNegWeight.assign(fEngine.fLogicVarToDefaultNegWeight.begin(),
                                  fEngine.fLogicVarToDefaultNegWeight.end());
    fVarToCurrentPosWeight.assign(fEngine.fLogicVarToDefaultPosWeight.begin(),
                                  fEngine.fLogicVarToDefaultPosWeight.end());
}

Evidence::Evidence (OnlineEngine& engine) : fEngine(engine){
    fVarToCurrentNegWeight.resize(
        fEngine.fLogicVarToDefaultNegWeight.size());
    fVarToCurrentPosWeight.resize(
        fEngine.fLogicVarToDefaultPosWeight.size());
    retractAll ();
}

void Evidence::varCommit (Variable v, int u) {
    varSet (v, 0.0);
    setCurrentWeightToDefault (
        fEngine.fSrcVarToSrcValToIndicator[v]->at(u));
}

void Evidence::valCommit (Variable v, int u, bool allow) {
    int l = fEngine.fSrcVarToSrcValToIndicator[v]->at(u);
    if (allow) {
        setCurrentWeightToDefault (l);
    } else {
        setCurrentWeight (l, 0.0);
    }
}

void Evidence::varRetract (Variable v) {
    setCurrentWeightsToDefaults (*(fEngine.fSrcVarToSrcValToIndicator[v]));
}

void Evidence::varSet (Variable v, double w) {
    setCurrentWeights (w, *(fEngine.fSrcVarToSrcValToIndicator[v]));
}

void Evidence::parmCommit (Potential t, int p, double w) {
    int l = fEngine.fSrcPotToSrcPosToParameter[t]->at(p);
    if (p == 0) {
        throw invalid_argument("Attempt to set value of parameter illegally!");
    }
    setCurrentWeight (l, w);
}

void Evidence::parmRetract (Potential t, int p) {
    int l = fEngine.fSrcPotToSrcPosToParameter[t]->at(p);
    if (p == 0) {
        throw invalid_argument("Attempt to set value of parameter illegally!");
    }
    setCurrentWeightToDefault (l);
}

vector<double> OnlineEngine::clone(vector<double> a) {
    vector<double> ans;
    ans.reserve(a.size());
    ans.assign(a.begin(), a.end());
    return ans;
}

void OnlineEngine::upwardPass (Evidence ev) {
    vector<double>& negValues = ev.fVarToCurrentNegWeight;
    vector<double>& posValues = ev.fVarToCurrentPosWeight;
    int numNodes = numAcNodes ();
    for (int n = 0; n < numNodes; n++) {
        char type = fNodeToType[n];
        if (type == MULTIPLY) {
            double v = 1.0;
            int last = fNodeToLastEdge[n];
            for (int e = first (n); e < last; e++) {
                int ch = fEdgeToTailNode[e];
                double chVal = fNodeToValue[ch];
                if (chVal == 0.0) {
                    v = 0.0;
                    break;
                }
                v *= chVal;
                if (v == 0.0) {
		  throw runtime_error("Underflow");
                }
            }
            fNodeToValue[n] = v;
        } else if (type == ADD) {
            double v = 0.0;
            int last = fNodeToLastEdge[n];
            for (int e = first (n); e < last; e++) {
                int ch = fEdgeToTailNode[e];
                v += fNodeToValue[ch];
            }
            fNodeToValue[n] = v;
        } else if (type == LITERAL) {
            int l = fNodeToLit[n];
            fNodeToValue[n] = (l < 0 ? negValues[-l] : posValues[l]);
        }
        // do nothing for a constant
    }
}

void OnlineEngine::twoPasses (Evidence ev) {

    // Upward pass.

    vector<double>& negValues = ev.fVarToCurrentNegWeight;
    vector<double>& posValues = ev.fVarToCurrentPosWeight;
    int numNodes = numAcNodes ();
    for (int n = 0; n < numNodes; n++) {
        fNodeToDerivative[n] = 0.0;
        char type = fNodeToType[n];
        if (type == MULTIPLY) {
            int numZeros = 0;
            double v = 1.0;
            int last = fNodeToLastEdge[n];
            for (int e = first (n); e < last; e++) {
                int ch = fEdgeToTailNode[e];
                double chVal = computedValue (ch);
                if (chVal == 0.0) {
                    if (++numZeros > 1) {
                        v = 0;
                        break;
                    }
                } else {
                    v *= chVal;
                    if (v == 0.0) {
                        throw runtime_error("Underflow");
                    }
                }
            }
            fNodeToValue[n] = v;
            fNodeToOneZero[n] = numZeros == 1;
        } else if (type == ADD) {
            double v = 0.0;
            int last = fNodeToLastEdge[n];
            for (int e = first (n); e < last; e++) {
                int ch = fEdgeToTailNode[e];
                double chVal = computedValue (ch);
                v += chVal;
            }
            fNodeToValue[n] = v;
        } else if (type == LITERAL) {
            int l = fNodeToLit[n];
            fNodeToValue[n] = (l < 0 ? negValues[-l] : posValues[l]);
        }
        // do nothing for a constant
    }

    // Downward pass.

    fNodeToDerivative[numNodes - 1] = 1.0;
    for (int n = numNodes - 1; n >= 0; n--) {
        char type = fNodeToType[n];
        if (type == LITERAL || type == CONSTANT) {
            continue;
        }
        int last = fNodeToLastEdge[n];
        if (type == MULTIPLY) {
            double value = fNodeToValue[n];
            if (value == 0.0) {
                continue;   // more than one zero
            }
            double x = fNodeToDerivative[n];
            if (x == 0.0) {
                continue;
            }
            x *= value;
            if (x == 0.0) {
                throw runtime_error("Underflow");
            }
            if (fNodeToOneZero[n]) { // exactly one zero
                for (int e = first (n); e < last; e++) {
                    int ch = fEdgeToTailNode[e];
                    double chVal = computedValue (ch);
                    if (chVal == 0.0) {
                        fNodeToDerivative[ch] += x;
                        break;
                    }
                }
            } else { // no zeros
                for (int e = first (n); e < last; e++) {
                    int ch = fEdgeToTailNode[e];
                    double chVal = computedValue (ch);
                    fNodeToDerivative[ch] += x / chVal;
                }
            }
        } else { /* PLUS NODE */
            double x = fNodeToDerivative[n];
            for (int e = first (n); e < last; e++) {
                int ch = fEdgeToTailNode[e];
                fNodeToDerivative[ch] += x;
            }
        }
    }

}

double OnlineEngine::computedValue (int n) {
    return fNodeToOneZero[n] ? 0 : fNodeToValue[n];
}

int OnlineEngine::first (int n) {
    return (n == 0) ? 0 : fNodeToLastEdge[n-1];
}

int OnlineEngine::rootNode () {
    return fNodeToType.size() - 1;
}

int OnlineEngine::numAcNodes () {
    return fNodeToType.size();
}


OnlineEngine::OnlineEngine (string acFilename, string lmFilename) {
	ifstream ac_fs (acFilename, ifstream::in);
        ifstream lm_fs (lmFilename, ifstream::in);
        initialize (ac_fs, lm_fs);
        ac_fs.close();
        lm_fs.close();
}

void OnlineEngine::initialize (istream& acReader, istream& lmReader) {
    readArithmeticCircuit(acReader);
    readLiteralMap(lmReader);
    fNodeToValue.resize(numAcNodes());
    fNodeToDerivative.resize(numAcNodes());
    fNodeToOneZero.resize(numAcNodes());
    fUpwardPassCompleted = false;
    fTwoPassesCompleted = false;
}

Variable OnlineEngine::varForName (string n) {
    return fNameToSrcVar[n];
}

Potential OnlineEngine::potForName (string n) {
    return fNameToSrcPot[n];
}

set<Variable> OnlineEngine::variables () {
    return fVariables;
}

set<Potential> OnlineEngine::potentials () {
    return fPotentials;
}

void OnlineEngine::assertEvidence (Evidence e, bool secondPass) {
    if (secondPass) {
        fAcVarToMostRecentNegWeight = clone (e.fVarToCurrentNegWeight);
        fAcVarToMostRecentPosWeight = clone (e.fVarToCurrentPosWeight);
        twoPasses (e);
    } else {
        fAcVarToMostRecentNegWeight.clear();
        fAcVarToMostRecentPosWeight.clear();
        upwardPass (e);
    }
    fUpwardPassCompleted = true;
    fTwoPassesCompleted = secondPass;
}

double OnlineEngine::probOfEvidence () {
    if (!fUpwardPassCompleted) {
        throw runtime_error ("assertEvidence () must be called!");
    }
    int root = rootNode ();
    return fTwoPassesCompleted ? computedValue (root) : fNodeToValue[root];
}

vector<double> OnlineEngine::varPartials (Variable v) {
    if (!fTwoPassesCompleted) {
        throw runtime_error (
            "assertEvidence () must be called with second pass flag set!");
    }
    vector<double> ans(v.domainNames().size());
    vector<int>& inds = *(fSrcVarToSrcValToIndicator[v]);
    for (int u = 0; u < ans.size(); u++) {
        int l = inds[u];
        ans[u] =
            (l < 0) ?
            fNodeToDerivative[fVarToNegLitNode[-l]] :
            fNodeToDerivative[fVarToPosLitNode[l]];
    }
    return ans;
}

map<Variable,vector<double> > OnlineEngine::varPartials (set<Variable> vs) {
    map<Variable,vector<double> > ans;
    for (Variable v : vs) {
        ans[v] = varPartials (v);
    }
    return ans;
}

vector<double> OnlineEngine::varMarginals (Variable v) {
    if (!fTwoPassesCompleted) {
        throw runtime_error (
            "assertEvidence () must be called with marginals flag set!");
    }
    vector<double> ans(v.domainNames ().size ());
    vector<int>& inds = *(fSrcVarToSrcValToIndicator[v]);
    for (int u = 0; u < ans.size(); u++) {
        int l = inds[u];
        ans[u] =
            (l < 0) ?
            fAcVarToMostRecentNegWeight[-l] *
            fNodeToDerivative[fVarToNegLitNode[-l]] :
            fAcVarToMostRecentPosWeight[l] *
            fNodeToDerivative[fVarToPosLitNode[l]];
    }
    return ans;
}

map<Variable,vector<double> > OnlineEngine::varMarginals (set<Variable> vs) {
    map<Variable,vector<double> > ans;
    for (Variable v : vs) {
        ans[v] = varMarginals (v);
    }
    return ans;
}

vector<double> OnlineEngine::varPosteriors (Variable v) {
    double pe = probOfEvidence ();
    vector<double> ans = varMarginals (v);
    for (int i = 0; i < ans.size(); i++) {
        ans[i] /= pe;
    }
    return ans;
}

map<Variable,vector<double> > OnlineEngine::varPosteriors (set<Variable> vs) {
    map<Variable,vector<double> > ans;
    for (Variable v : vs) {
        ans[v] = varPosteriors (v);
    }
    return ans;
}

vector<double> OnlineEngine::potPartials (Potential pot) {
    vector<int>& parms = *(fSrcPotToSrcPosToParameter[pot]);
    vector<double> ans(parms.size());
    for (int pos = 0; pos < ans.size(); pos++) {
        int l = parms[pos];
        ans[pos] =
            l == 0  ? NAN :
            l < 0 ? fNodeToDerivative[fVarToNegLitNode[-l]] :
            fNodeToDerivative[fVarToPosLitNode[l]];
    }
    return ans;
}

map<Potential,vector<double> > OnlineEngine::potPartials (
    set<Potential> ps) {
    map<Potential,vector<double> > ans;
    for (Potential p : ps) {
        ans[p] = potPartials (p);
    }
    return ans;
}

vector<double> OnlineEngine::potMarginals (Potential p) {
    vector<int>& parms = *(fSrcPotToSrcPosToParameter[p]);
    vector<double> ans = potPartials (p);
    for (int pos = 0; pos < ans.size(); pos++) {
        int l = parms[pos];
        if (!isnan (ans[pos])) {
            ans[pos] *=
                l < 0 ?
                fAcVarToMostRecentNegWeight[-l] : fAcVarToMostRecentPosWeight[l];
        }
    }
    return ans;
}

map<Potential,vector<double> > OnlineEngine::potMarginals (set<Potential> ps) {
    map<Potential,vector<double> > ans;
    for (Potential p : ps) {
        ans[p] = potMarginals (p);
    }
    return ans;
}

vector<double> OnlineEngine::potPosteriors (Potential p) {
    vector<double> ans = potMarginals (p);
    double pe = probOfEvidence ();
    for (int pos = 0; pos < ans.size(); pos++) {
        if (!isnan (ans[pos])) {
            ans[pos] /= pe;
        }
    }
    return ans;
}

map<Potential,vector<double> > OnlineEngine::potPosteriors (
    set<Potential> ps) {
    map<Potential,vector<double> > ans;
    for (Potential p : ps) {
        ans[p] = potPosteriors (p);
    }
    return ans;
}

const string OnlineEngine::READ_DELIMITER  = "\\$";
const string OnlineEngine::DELIMITER = "$";

void OnlineEngine::readArithmeticCircuit(istream& r) {

    int numNodes = __INT_MAX__;
    int nextEdge = 0;
    int nextNode = 0;
    
    // Process each line.
    
    while (nextNode < numNodes) {
      
      // Read the line.  Quit if eof.  Skip if comment or blank line.
      // Otherwise, split into tokens.
      
      string line;
      if (!getline(r, line)) {break;} // eof
      if (boost::starts_with(line, "c")) {continue;} // comment
      boost::trim(line);
      if (line.length () == 0) {continue;} // blank line
      vector<string> tokens;
      boost::split(tokens, line, boost::is_any_of("\t "), boost::token_compress_on);
      
      // A header line looks like: "nnf" numNodes numEdges numVars
      
      if (tokens[0] == "nnf") {
	numNodes = boost::lexical_cast<int>(tokens[1]);
        int numEdges = boost::lexical_cast<int>(tokens[2]);
        int numAcVars = boost::lexical_cast<int>(tokens[3]);
        fEdgeToTailNode.resize(numEdges);
        fNodeToLastEdge.resize(numNodes);
        fNodeToType.resize(numNodes);
        fNodeToLit.resize(numNodes, 0);
        fVarToNegLitNode.assign(numAcVars + 1, -1);
        fVarToPosLitNode.assign(numAcVars + 1, -1);
        continue;
      }
      
      // This is not a header line, so it must be a node line, which looks
      // like one of the following:
      //   "A" numChildren child+
      //   "O" splitVar numChildren child+
      //   "L" literal
      char ch = tokens[0].at(0);
      if (ch == 'A') {
        fNodeToType[nextNode] = MULTIPLY;
        for (int chIndex = 2; chIndex < tokens.size(); chIndex++) {
          fEdgeToTailNode[nextEdge++] = boost::lexical_cast<int> (tokens[chIndex]);
        }
      } else if (ch == 'O') {
        fNodeToType[nextNode] = ADD;
        for (int chIndex = 3; chIndex < tokens.size(); chIndex++) {
          fEdgeToTailNode[nextEdge++] = boost::lexical_cast<int> (tokens[chIndex]);
        }
      } else /* ch == 'L' */ {
        fNodeToType[nextNode] = LITERAL;
        int l = boost::lexical_cast<int> (tokens[1]);
        fNodeToLit[nextNode] = l;
        (l < 0 ? fVarToNegLitNode : fVarToPosLitNode)[abs (l)] =
          nextNode;
      }
      fNodeToLastEdge[nextNode] = nextEdge;
      nextNode++;

    }

}

void OnlineEngine::readLiteralMap(istream& r) {
    
    // Prepare to parse.
    
    int numLits = __INT_MAX__;
    int litsFinished = 0;

    // Process each line.
    
    while (litsFinished < numLits) {
      
      // Read the line.  Quit if eof.  Skip if comment (including blank
      // lines).  Otherwise, split into tokens.
      
      
      string line;
      if (!getline(r, line)) {break;} // eof
      if (!boost::starts_with(line, ("cc" + DELIMITER))) {continue;} // comment
      boost::trim(line);
      vector<string> tokens;
      boost::split(tokens, line, boost::is_any_of(READ_DELIMITER), boost::token_compress_on);
      
      // If the line is a header, it is of the form: "cc" "N" numLogicVars.
      string type = tokens[1];
      if (type == "N") {
	int n = boost::lexical_cast<int>(tokens[2]);
        numLits = n * 2;
        fLogicVarToDefaultNegWeight.resize(n+1);
        fLogicVarToDefaultPosWeight.resize(n+1);
        continue;
      }
      
      // If the line is a variable line, then it is of the form:
      // "cc" "V" srcVarName numSrcVals (srcVal)+

      if (type == "V") {
        string vn = tokens[2];
        int valCount = boost::lexical_cast<int>(tokens[3]);
        vector<string> valNames(valCount, "");
        for (int i = 0; i < valCount; i++) {valNames[i] = tokens[4 + i];}
        Variable v(vn, valNames);
        fSrcVarToSrcValToIndicator[v] = new vector<int>(valCount);
        fNameToSrcVar[vn] = v;
        continue;
      }
      
      // If the line is a potential line, then it is of the form:
      // "cc" "T" srcPotName parameterCnt.
      
      if (type == "T") {
        string tn = tokens[2];
        int parmCount = boost::lexical_cast<int> (tokens[3]);
        Potential pot(tn, parmCount);
        fSrcPotToSrcPosToParameter[pot] = new vector<int>(parmCount);
        fNameToSrcPot[tn] = pot;
        continue;
      }

      // This is not a header line, a variable line, or potential line, so
      // it must be a literal description, which looks like one of the
      // following:
      //   "cc" "I" literal weight srcVarName srcValName srcVal
      //   "cc" "P" literal weight srcPotName pos+
      //   "cc" "A" literal weight
      
      int l = boost::lexical_cast<int>(tokens[2]);
      double w = boost::lexical_cast<double>(tokens[3]);
      (l < 0 ? fLogicVarToDefaultNegWeight : fLogicVarToDefaultPosWeight)
        [abs (l)] = w;
      if (type == "I") {
        string vn = tokens[4];
        //String un = tokens[5];
        int u = boost::lexical_cast<int>(tokens[6]);
        fSrcVarToSrcValToIndicator[fNameToSrcVar[vn]]->at(u) = l;
      } else if (type == "P") {
        string tn = tokens[4];
	vector<string> posStrings;
	boost::split(posStrings, tokens[5], boost::is_any_of(","));
        if (posStrings.size() == 1) {
	  int pos = boost::lexical_cast<int>(posStrings[0]);
          fSrcPotToSrcPosToParameter[fNameToSrcPot[tn]]->at(pos) = l;
        }
      } else if (type == "A") {
      } else {
        throw runtime_error (
          "\"cc\" must be followed by \"N\", \"V\", \"T\", \"I\", \"P\", or \"A\"");
      }
      ++litsFinished;
      
    }

    // Now create the variables, the map from variable name to variable, and
    // the map from variable to value to indicator.
    for (auto p : fNameToSrcVar)
      fVariables.insert(p.second);
    for (auto p : fNameToSrcPot)
      fPotentials.insert(p.second);

}
