#ifndef OCCA_PARSER_DEFINES_HEADER
#define OCCA_PARSER_DEFINES_HEADER

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <stack>
#include <map>

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#if (OCCA_OS & (LINUX_OS | OSX_OS))
#  include <unistd.h>
#else
#  include <io.h>
#endif

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "occaDefines.hpp"

namespace occa {
  namespace parserNS {
    class opHolder;
    class typeHolder;

    class macroInfo;

    class expNode;
    class statement;

    class attribute_t;
    class overloadedOp_t;

    class typeInfo;
    class varInfo;
    class kernelInfo;

    class parserBase;

    template <class TM> class node;
  };

  //---[ Info ]-----------------------------------
  typedef uint64_t info_t;

  typedef std::vector<int>                                  intVector_t;
  typedef std::vector<intVector_t>                          intVecVector_t;
  typedef std::vector<std::string>                          stringVector_t;
  typedef parserNS::node<parserNS::statement*>              statementNode;
  typedef parserNS::node<parserNS::varInfo*>                varInfoNode;

  typedef std::vector<parserNS::statement*>                 statementVector_t;
  typedef std::vector<statementNode*>                       statementNodeVector_t;
  typedef std::vector<parserNS::varInfo*>                   varInfoVector_t;
  typedef std::vector<varInfoVector_t>                      varInfoVecVector_t;

  typedef std::map<std::string,info_t>                      macroMap_t;
  typedef macroMap_t::iterator                              macroMapIterator;
  typedef macroMap_t::const_iterator                        cMacroMapIterator;

  typedef std::map<std::string,info_t>                      keywordTypeMap_t;
  typedef keywordTypeMap_t::iterator                        keywordTypeMapIterator;
  typedef keywordTypeMap_t::const_iterator                  cKeywordTypeMapIterator;

  typedef std::map<parserNS::opHolder,info_t>               opTypeMap_t;
  typedef opTypeMap_t::iterator                             opTypeMapIterator;
  typedef opTypeMap_t::const_iterator                       cOpTypeMapIterator;

  typedef std::vector<parserNS::typeInfo*>                  anonymousTypeMap_t;

  typedef std::map<std::string,parserNS::typeInfo*>         scopeTypeMap_t;
  typedef scopeTypeMap_t::iterator                          scopeTypeMapIterator;
  typedef scopeTypeMap_t::const_iterator                    cScopeTypeMapIterator;

  typedef std::map<std::string,parserNS::varInfo*>          scopeVarMap_t;
  typedef scopeVarMap_t::iterator                           scopeVarMapIterator;
  typedef scopeVarMap_t::const_iterator                     cScopeVarMapIterator;

  typedef std::map<std::string,parserNS::overloadedOp_t*>   opOverloadMaps_t;
  typedef opOverloadMaps_t::iterator                        opOverloadMapsIterator;
  typedef opOverloadMaps_t::const_iterator                  cOpOverloadMapsIterator;

  typedef std::map<parserNS::varInfo*,parserNS::statement*> varOriginMap_t;
  typedef varOriginMap_t::iterator                          varOriginMapIterator;
  typedef varOriginMap_t::const_iterator                    cVarOriginMapIterator;

  typedef std::vector<varOriginMap_t>                       varOriginMapVector_t;
  typedef varOriginMapVector_t::iterator                    varOriginMapVectorIterator;
  typedef varOriginMapVector_t::const_iterator              cVarOriginMapVectorIterator;

  typedef std::map<parserNS::varInfo*, statementNode>       varUsedMap_t;
  typedef varUsedMap_t::iterator                            varUsedMapIterator;
  typedef varUsedMap_t::const_iterator                      cVarUsedMapIterator;

  typedef std::vector<varUsedMap_t>                         varUsedMapVector_t;
  typedef varUsedMapVector_t::iterator                      varUsedMapVectorIterator;
  typedef varUsedMapVector_t::const_iterator                cVarUsedMapVectorIterator;

  typedef std::map<std::string,parserNS::kernelInfo*>       kernelInfoMap_t;
  typedef kernelInfoMap_t::iterator                         kernelInfoIterator;
  typedef kernelInfoMap_t::const_iterator                   cKernelInfoIterator;

  typedef std::map<std::string, std::string>                strToStrMap_t;
  typedef strToStrMap_t::iterator                           strToStrMapIterator;
  typedef strToStrMap_t::const_iterator                     cStrToStrMapIterator;

  typedef std::map<std::string,parserNS::attribute_t*>      attributeMap_t;
  typedef attributeMap_t::iterator                          attributeMapIterator;
  typedef attributeMap_t::const_iterator                    cAttributeMapIterator;

  typedef std::map<parserNS::varInfo*,parserNS::varInfo*>   varToVarMap_t;
  typedef varToVarMap_t::iterator                           varToVarMapIterator;
  typedef varToVarMap_t::const_iterator                     cVarToVarMapIterator;

  typedef std::map<parserNS::statement*,int>                statementIdMap_t;
  typedef statementIdMap_t::iterator                        statementIdMapIterator;

  typedef std::map<parserNS::varInfo*,int>                  varInfoIdMap_t;
  typedef varInfoIdMap_t::iterator                          varInfoIdMapIterator;

  typedef std::map<int,bool>                                idDepMap_t;
  typedef idDepMap_t::iterator                              idDepMapIterator;

  typedef std::vector<parserNS::expNode*>                   expVector_t;
  typedef expVector_t::iterator                             expVectorIterator;

  typedef void (parserNS::parserBase::*applyToAllStatements_t)(parserNS::statement &s);

  namespace parserInfo {
    static const int nothing            = 0;
    static const int parsingC           = (1 << 0);
    static const int parsingFortran     = (1 << 1);

    //---[ Check Flags ]------
    static const int checkSubStatements = (1 << 0);
  };

  namespace parserNS {

    extern keywordTypeMap_t *keywordType;
    extern keywordTypeMap_t cKeywordType, fortranKeywordType;

    extern bool cKeywordsAreInitialized;
    extern bool fortranKeywordsAreInitialized;

    //   ---[ Delimiters ]---------
    static const char whitespace[]      = " \t\r\n\v\f\0";

    static const char cWordDelimiter[]  = " \t\r\n\v\f!\"#%&'()*+,-./:;<=>?[]^{|}~@\0";
    static const char cWordDelimiter2[] = "!=###<#>%=&&&=*=+++=-=--->.*../=::<<<===>=>>^=|=||\0";
    static const char cWordDelimiter3[] = "->*...<<=>>=\0";

    static const char fortranWordDelimiter[]  = " \t\r\n\v\f\"#%'()*+,-./;<=>\0";
    static const char fortranWordDelimiter2[] = "**/=::<===>=\0";

    //   ---[ Types ]-----------
    static const info_t isTypeDef = (1 << 0);
    static const info_t isVarInfo = (1 << 1);

    static const info_t structType          = (1 <<  0);
    static const info_t noType              = (1 <<  1);
    static const info_t boolType            = (1 <<  2);
    static const info_t charType            = (1 <<  3);
    static const info_t ushortType          = (1 <<  4);
    static const info_t shortType           = (1 <<  5);
    static const info_t uintType            = (1 <<  6);
    static const info_t intType             = (1 <<  7);
    static const info_t ulongType           = (1 <<  8);
    static const info_t longType            = (1 <<  9);
    static const info_t ulonglongType       = (1 << 10);
    static const info_t longlongType        = (1 << 11);
    static const info_t voidType            = (1 << 12);
    static const info_t floatType           = (1 << 13);
    static const info_t doubleType          = (1 << 14);

    //   ---[ Instructions ]----
    static const info_t doNothing           = (1 << 0);
    static const info_t ignoring            = (3 << 1);
    static const info_t ignoreUntilNextHash = (1 << 1);
    static const info_t ignoreUntilEnd      = (1 << 2);
    static const info_t readUntilNextHash   = (1 << 3);
    static const info_t doneIgnoring        = (1 << 4);
    static const info_t startHash           = (1 << 5);

    static const info_t keepMacro           = (1 << 6);
    static const info_t forceLineRemoval    = (1 << 7);

    static const info_t readingCode          = 0;
    static const info_t insideCommentBlock   = 1;
    static const info_t finishedCommentBlock = 2;
    //==============================================
  };
};

#endif
