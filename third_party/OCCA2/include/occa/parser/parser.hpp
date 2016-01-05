#ifndef OCCA_PARSER_HEADER
#define OCCA_PARSER_HEADER

#include "occa/parser/defines.hpp"
#include "occa/parser/preprocessor.hpp"
#include "occa/parser/tools.hpp"
#include "occa/parser/nodes.hpp"
#include "occa/parser/types.hpp"
#include "occa/parser/statement.hpp"
#include "occa/parser/analyzer.hpp"
#include "occa/parser/magic.hpp"
#include "occa/tools.hpp"

namespace occa {
  namespace parserNS {
    class occaLoopInfo;

    extern intVector_t loadedLanguageVec;

    int loadedLanguage();

    void pushLanguage(const int language);
    int popLanguage();

    class parserBase {
    public:
      std::string filename;

      int parsingLanguage;

      flags_t parsingFlags;

      macroMap_t macroMap;
      std::vector<macroInfo> macros;

      //---[ Parser Warnings ]----------
      bool macrosAreInitialized;
      bool _hasMagicEnabled;
      bool _compilingForCPU;
      bool _warnForMissingBarriers;
      bool _warnForConditionalBarriers;
      bool _insertBarriersAutomatically;
      //================================

      varOriginMap_t varOriginMap;

      kernelInfoMap_t kernelInfoMap;

      statement *globalScope;

      parserBase();
      ~parserBase();

      inline const std::string parseFile(const std::string &filename_,
                                         const flags_t &flags_ = flags_t()){

        return parseFile("", filename_, flags_);
      }

      const std::string parseFile(const std::string &header,
                                  const std::string &filename,
                                  const flags_t &flags_ = flags_t());

      const std::string parseSource(const char *cRoot);

      //---[ Parser Warnings ]----------
      void loadParserFlags(const flags_t &flags);

      bool hasMagicEnabled();
      bool compilingForCPU();
      bool warnForMissingBarriers();
      bool warnForConditionalBarriers();
      bool insertBarriersAutomatically();
      //================================

      //---[ Macro Parser Functions ]---
      std::string getMacroName(const char *&c);
      std::string getMacroIncludeFile(const char *&c);

      typeHolder evaluateMacroStatement(const char *c);
      bool evaluateMacroBoolStatement(const char *c);

      void loadMacroInfo(macroInfo &info, const char *&c);

      int loadMacro(expNode &allExp, int leafPos, const int state = doNothing);
      int loadMacro(const std::string &line, const int state = doNothing);
      int loadMacro(expNode &allExp, int leafPos, const std::string &line, const int state = doNothing);

      void applyMacros(std::string &line);

      void preprocessMacros(expNode &allExp);

      expNode splitAndPreprocessContent(const std::string &s,
                                        const int parsingLanguage_ = parserInfo::parsingC);
      expNode splitAndPreprocessContent(const char *cRoot,
                                        const int parsingLanguage_ = parserInfo::parsingC);
      //====================================

      void initMacros(const int parsingLanguage_ = parserInfo::parsingC);
      void initFortranMacros();

      void loadLanguageTypes();

      void applyToAllStatements(statement &s,
                                applyToAllStatements_t func);

      void applyToAllKernels(statement &s,
                             applyToAllStatements_t func);

      static bool statementIsAKernel(statement &s);

      static statement* getStatementKernel(statement &s);

      bool statementKernelUsesNativeOCCA(statement &s);

      bool statementKernelUsesNativeOKL(statement &s);

      bool statementKernelUsesNativeLanguage(statement &s);


      void updateOccaOMLoopAttributes(statement &s,
                                      const std::string &loopTag,
                                      const std::string &loopNest);

      void setupOccaFors(statement &s);

      void removeIntraLoopDepsFromIterExp(statement &s);

      bool statementIsOccaOuterFor(statement &s);
      bool statementIsOccaInnerFor(statement &s);

      bool statementHasOccaOuterFor(statement &s);
      bool statementHasOccaFor(statement &s);

      bool statementHasOklFor(statement &s);

      bool statementHasOccaStuff(statement &s);

      //   ---[ Loop Reordering ]-------
      void reorderLoops();

      void reorderLoops(statementVector_t &loopsToReorder,
                        const int start,
                        const int end);

      intVector_t relatedReorderLoops(statementVector_t &loopsToReorder,
                                      const int start,
                                      const int end);

      void placeLoopsToReorder(statement &s,
                               statementVector_t &loopsToReorder);
      //   =============================

      void retagOccaLoops();
      void retagOccaLoops(statement &s);

      void splitTileOccaFors(statement &s);

      void markKernelFunctions();

      void labelNativeKernels();

      void setupCudaVariables(statement &s);

      void addFunctionPrototypes();

      void checkOccaBarriers(statement &s);
      void addOccaBarriers();
      void addOccaBarriersToStatement(statement &s);

      bool statementHasBarrier(statement &s);

      void addParallelFors(statement &s);

      void updateConstToConstant();

      void addArgQualifiers();
      void addArgQualifiersTo(statement &s);

      void floatSharedAndExclusivesUp(statement &s);
      statementNode* appendSharedAndExclusives(statement &s,
                                               statementNode *snTail,
                                               bool isAppending = false);

      void modifyExclusiveVariables(statement &s);

      void modifyTextureVariables();

      //   ---[ Load Kernels ]----------
      void loadKernelInfos();

      statementNode* splitKernelStatement(statementNode *snKernel,
                                          kernelInfo &info);

      statementVector_t findOuterLoopSets(statement &sKernel);
      void findOuterLoopSets(statement &s, statementVector_t &omLoops);

      statementVector_t findOccaLoops(statement &sKernel);
      void findOccaLoops(statement &s, statementVector_t &occaLoops);

      varOriginMap_t findKernelDependenciesFor(statement &sKernel,
                                               statement &omLoop);

      varOriginMap_t findDependenciesFor(statement &s,
                                         const int flags = parserInfo::checkSubStatements);

      void findDependenciesFor(statement &s,
                               varOriginMap_t &deps,
                               const int flags = parserInfo::checkSubStatements);

      varOriginMap_t findDependenciesFor(expNode &e);

      void findDependenciesFor(expNode &e,
                               varOriginMap_t &deps);

      statementVector_t newKernelsFromLoops(statement &sKernel,
                                            statementVector_t &omLoops,
                                            varOriginMapVector_t &varDeps);

      void addDepStatementsToKernel(statement &sKernel,
                                    varOriginMap_t &deps);

      void addDepsToKernelArguments(statement &sKernel,
                                    varOriginMap_t &deps);

      statement& launchStatementForKernel(statement &sKernel,
                                          statement &omLoop,
                                          const int newKernelPos,
                                          varInfo &newKernelVar);

      void storeKernelInfo(kernelInfo &info,
                           statement &sKernel,
                           statementVector_t &newKernels);

      void zeroOccaIdsFrom(statement &s);
      void zeroOccaIdsFrom(expNode &e);

      void addNestedKernelArgTo(statement &sKernel);
      //   =============================

      static int getKernelOuterDim(statement &s);
      static int getKernelInnerDim(statement &s);
      static int getKernelDimFor(statement &s, const std::string &tag);

      int getOuterMostForDim(statement &s);
      int getInnerMostForDim(statement &s);
      int getForDim(statement &s, const std::string &tag);

      void splitDefineForVariable(varInfo &var);
      void splitDefineAndInitForVariable(varInfo &var);

      void addInnerFors(statement &s);
      void addInnerForsTo(statement &s,
                          varInfoIdMap_t &varInfoIdMap,
                          int &currentInnerID,
                          const int innerDim);

      void checkStatementForExclusives(statement &s,
                                       varInfoIdMap_t &varInfoIdMap,
                                       const int innerID);

      void addOuterFors(statement &s);

      void removeUnnecessaryBlocksInKernel(statement &s);
      void addOccaForsToKernel(statement &s);

      void addOccaFors();

      void setupOccaVariables(statement &s);

      //---[ Operator Information ]---------------
      varInfo* hasOperator(const info_t expInfo,
                           const std::string &op,
                           varInfo &var);

      varInfo* hasOperator(const info_t expInfo,
                           const std::string &op,
                           varInfo &varL,
                           varInfo &varR);

      varInfo thVarInfo(const info_t thType);

      varInfo thOperatorReturnType(const info_t expInfo,
                                   const std::string &op,
                                   const info_t thType);

      varInfo thOperatorReturnType(const info_t expInfo,
                                   const std::string &op,
                                   const info_t thTypeL,
                                   const info_t thTypeR);
      //==========================================
    };

    bool isAnOccaID(const std::string &s);
    bool isAnOccaInnerID(const std::string &s);
    bool isAnOccaOuterID(const std::string &s);
    bool isAnOccaGlobalID(const std::string &s);

    bool isAnOccaDim(const std::string &s);
    bool isAnOccaInnerDim(const std::string &s);
    bool isAnOccaOuterDim(const std::string &s);
    bool isAnOccaGlobalDim(const std::string &s);

    expNode splitContent(const std::string &str,
                         const int parsingLanguage_ = parserInfo::parsingC);
    expNode splitContent(const char *cRoot,
                         const int parsingLanguage_ = parserInfo::parsingC);

    expNode splitAndLabelContent(const std::string &str,
                                 const int parsingLanguage_ = parserInfo::parsingC);
    expNode splitAndLabelContent(const char *cRoot,
                                 const int parsingLanguage_ = parserInfo::parsingC);
    expNode splitAndOrganizeContent(const std::string &str,
                                    const int parsingLanguage_ = parserInfo::parsingC);
    expNode splitAndOrganizeContent(const char *cRoot,
                                    const int parsingLanguage_ = parserInfo::parsingC);

    expNode& labelCode(expNode &allExp,
                       const int parsingLanguage_ = parserInfo::parsingC);

    bool checkLastTwoNodes(expNode &node,
                           const std::string &leftValue,
                           const std::string &rightValue,
                           const int parsingLanguage_ = parserInfo::parsingC);

    void mergeLastTwoNodes(expNode &node,
                           const bool addSpace = true,
                           const int parsingLanguage_ = parserInfo::parsingC);

    expNode createExpNodeFrom(const std::string &source);
    expNode createOrganizedExpNodeFrom(const std::string &source);

    void loadKeywords(const int parsingLanguage_ = parserInfo::parsingC);
    void loadCKeywords();
    void loadFortranKeywords();
    void initCKeywords();
    void initFortranKeywords();

    //---[ OCCA Loop Info ]-------------
    class occaLoopInfo {
    public:
      statement *sInfo;
      int parsingLanguage;

      occaLoopInfo(statement &s,
                   const int parsingLanguage_ = parserInfo::parsingC,
                   const std::string &tag = "");

      void lookForLoopFrom(statement &s,
                           const std::string &tag = "");

      void loadForLoopInfo(int &innerDims, int &outerDims,
                           std::string *innerIters,
                           std::string *outerIters);

      void getLoopInfo(std::string &loopTag,
                       std::string &loopNest);

      void getLoopNode1Info(std::string &iter,
                            std::string &start);

      void getLoopNode2Info(std::string &bound,
                            std::string &iterCheck);

      void getLoopNode3Info(std::string &stride,
                            std::string &strideOpSign,
                            std::string &strideOp);

      void setIterDefaultValues();

      std::string getSetupExpression();
    };
    //==================================
    //==============================================
  }

  // Just to ignore the namespace
  class parser: public parserNS::parserBase {};
}

#endif
