#include "occa.hpp"
#include "occaParser.hpp"

// Use events for timing!

namespace occa {
  //---[ Globals & Flags ]------------
  const int parserVersion = 100;

  kernelInfo defaultKernelInfo;

  const int autoDetect = (1 << 0);
  const int srcInUva   = (1 << 1);
  const int destInUva  = (1 << 2);

  bool uvaEnabledByDefault_f = false;
  bool verboseCompilation_f  = true;

  void setVerboseCompilation(const bool value){
    verboseCompilation_f = value;
  }
  //==================================


  //---[ UVA ]------------------------
  ptrRangeMap_t uvaMap;
  memoryVector_t uvaDirtyMemory;

  bool hasUvaEnabledByDefault(){
    return uvaEnabledByDefault_f;
  }

  void enableUvaByDefault(){
    uvaEnabledByDefault_f = true;
  }

  void disableUvaByDefault(){
    uvaEnabledByDefault_f = false;
  }

  ptrRange_t::ptrRange_t() :
    start(NULL),
    end(NULL) {}

  ptrRange_t::ptrRange_t(void *ptr, const uintptr_t bytes) :
    start((char*) ptr),
    end(((char*) ptr) + bytes) {}

  ptrRange_t::ptrRange_t(const ptrRange_t &r) :
    start(r.start),
    end(r.end) {}

  ptrRange_t& ptrRange_t::operator = (const ptrRange_t &r){
    start = r.start;
    end   = r.end;

    return *this;
  }

  bool ptrRange_t::operator == (const ptrRange_t &r) const {
    return ((start <= r.start) && (r.start < end));
  }

  bool ptrRange_t::operator != (const ptrRange_t &r) const {
    return ((r.start < start) || (end <= r.start));
  }

  int operator < (const ptrRange_t &a, const ptrRange_t &b){
    return ((a != b) && (a.start < b.start));
  }

  uvaPtrInfo_t::uvaPtrInfo_t() :
    mem(NULL) {}

  uvaPtrInfo_t::uvaPtrInfo_t(void *ptr){
    ptrRangeMap_t::iterator it = uvaMap.find(ptr);

    if(it != uvaMap.end())
      mem = (it->second);
    else
      mem = (occa::memory_v*) ptr; // Defaults to ptr being a memory_v
  }

  uvaPtrInfo_t::uvaPtrInfo_t(occa::memory_v *mem_) :
    mem(mem_) {}

  uvaPtrInfo_t::uvaPtrInfo_t(const uvaPtrInfo_t &upi) :
    mem(upi.mem) {}

  uvaPtrInfo_t& uvaPtrInfo_t::operator = (const uvaPtrInfo_t &upi){
    mem = upi.mem;

    return *this;
  }

  occa::device uvaPtrInfo_t::getDevice(){
    occa::memory m(mem);

    return occa::device(m.getOccaDeviceHandle());
  }

  occa::memory uvaPtrInfo_t::getMemory(){
    return occa::memory(mem);
  }

  occa::memory_v* uvaToMemory(void *ptr){
    ptrRangeMap_t::iterator it = uvaMap.find(ptr);

    if(it == uvaMap.end())
      return NULL;

    return it->second;
  }

  void syncToDevice(void *ptr, const uintptr_t bytes){
    occa::memory_v *mem = uvaToMemory(ptr);

    if(mem == NULL)
      return;

    if(!mem->dHandle->fakesUva())
      memcpy(mem->handle, mem->uvaPtr, bytes);
    else
      occa::memory(mem).syncToDevice(bytes);
  }

  void syncFromDevice(void *ptr, const uintptr_t bytes){
    occa::memory_v *mem = uvaToMemory(ptr);

    if(mem == NULL)
      return;

    if(!mem->dHandle->fakesUva())
      memcpy(mem->uvaPtr, mem->handle, bytes);
    else
      occa::memory(mem).syncFromDevice(bytes);
  }

  bool needsSync(void *ptr){
    occa::memory m(ptr);

    return m.uvaIsDirty();
  }

  void dontSync(void *ptr){
    removeFromDirtyMap(ptr);
  }

  void removeFromDirtyMap(void *ptr){
    ptrRangeMap_t::iterator it = uvaMap.find(ptr);

    if(it == uvaMap.end())
      return;

    memory m(it->second);

    if(!m.uvaIsDirty())
      return;

    removeFromDirtyMap(m.getOccaMemoryHandle());
  }

  void removeFromDirtyMap(memory_v *mem){
    occa::memory m(mem);

    const size_t dirtyEntries = uvaDirtyMemory.size();

    for(size_t i = 0; i < dirtyEntries; ++i){
      if(uvaDirtyMemory[i] == mem){
        m.uvaMarkClean();
        uvaDirtyMemory.erase(uvaDirtyMemory.begin() + i);

        break;
      }
    }
  }

  void setupMagicFor(void *ptr){
    ptrRangeMap_t::iterator it = uvaMap.find(ptr);

    if(it == uvaMap.end())
      return;

    memory_v &mem = *(it->second);

    if(mem.dHandle->fakesUva())
      return;

    if(mem.uvaPtr == NULL)
      mem.uvaPtr = cpu::malloc(mem.size);

    memcpy(mem.uvaPtr, mem.handle, mem.size);
  }

  void free(void *ptr){
    ptrRangeMap_t::iterator it = uvaMap.find(ptr);

    if((it != uvaMap.end()) &&
       (((void*) it->first.start) != ((void*) it->second))){

      occa::memory(it->second).free();
    }
    else
      ::free(ptr);
  }
  //==================================

  //---[ Helper Classes ]-------------
  const int uint8FormatIndex  = 0;
  const int uint16FormatIndex = 1;
  const int uint32FormatIndex = 2;
  const int int8FormatIndex   = 3;
  const int int16FormatIndex  = 4;
  const int int32FormatIndex  = 5;
  const int halfFormatIndex   = 6;
  const int floatFormatIndex  = 7;

  const int sizeOfFormats[8] = {1, 2, 4,
                                1, 2, 4,
                                2, 4};

  formatType::formatType(const int format__, const int count__){
    format_ = format__;
    count_  = count__;
  }

  formatType::formatType(const formatType &ft){
    format_ = ft.format_;
    count_  = ft.count_;
  }

  formatType& formatType::operator = (const formatType &ft){
    format_ = ft.format_;
    count_  = ft.count_;

    return *this;
  }

  int formatType::count() const {
    return count_;
  }

  size_t formatType::bytes() const {
    return (sizeOfFormats[format_] * count_);
  }

  const int readOnly  = 1;
  const int readWrite = 2;

  const occa::formatType uint8Format(uint8FormatIndex  , 1);
  const occa::formatType uint8x2Format(uint8FormatIndex, 2);
  const occa::formatType uint8x4Format(uint8FormatIndex, 4);

  const occa::formatType uint16Format(uint16FormatIndex  , 1);
  const occa::formatType uint16x2Format(uint16FormatIndex, 2);
  const occa::formatType uint16x4Format(uint16FormatIndex, 4);

  const occa::formatType uint32Format(uint32FormatIndex  , 1);
  const occa::formatType uint32x2Format(uint32FormatIndex, 2);
  const occa::formatType uint32x4Format(uint32FormatIndex, 4);

  const occa::formatType int8Format(int8FormatIndex  , 1);
  const occa::formatType int8x2Format(int8FormatIndex, 2);
  const occa::formatType int8x4Format(int8FormatIndex, 4);

  const occa::formatType int16Format(int16FormatIndex  , 1);
  const occa::formatType int16x2Format(int16FormatIndex, 2);
  const occa::formatType int16x4Format(int16FormatIndex, 4);

  const occa::formatType int32Format(int32FormatIndex  , 1);
  const occa::formatType int32x2Format(int32FormatIndex, 2);
  const occa::formatType int32x4Format(int32FormatIndex, 4);

  const occa::formatType halfFormat(halfFormatIndex  , 1);
  const occa::formatType halfx2Format(halfFormatIndex, 2);
  const occa::formatType halfx4Format(halfFormatIndex, 4);

  const occa::formatType floatFormat(floatFormatIndex  , 1);
  const occa::formatType floatx2Format(floatFormatIndex, 2);
  const occa::formatType floatx4Format(floatFormatIndex, 4);

  //---[ Arg Info Map ]-------
  argInfoMap::argInfoMap(){}

  argInfoMap::argInfoMap(const std::string &infos){
    if(infos.size() == 0)
      return;

    parserNS::strNode *n;

    n = parserNS::splitContent(infos);
    n = parserNS::labelCode(n);

    while(n){
      std::string &info = n->value;
      std::string value;

      n = n->right;

      if((info != "mode")        &&
         (info != "UVA")         &&
         (info != "platformID")  &&
         (info != "deviceID")    &&
         (info != "schedule")    &&
         (info != "chunk")       &&
         (info != "threadCount") &&
         (info != "schedule")    &&
         (info != "pinnedCores")){

        std::cout << "Flag [" << info << "] is not available, skipping it\n";

        while(n && (n->value != ","))
          n = n->right;

        if(n)
          n = n->right;

        continue;
      }

      if(n == NULL)
        break;

      if(n->value == "=")
        n = n->right;

      while(n && (n->value != ",")){
        std::string &v = n->value;

        occa::strip(v);

        if(v.size()){
          if(segmentPair(v[0]) == 0){
            value += v;
            value += ' ';
          }
          else if(n->down){
            std::string dv = n->down->toString();
            occa::strip(dv);

            value += dv;
            value += ' ';
          }
        }

        n = n->right;
      }

      if(n)
        n = n->right;

      occa::strip(value);

      iMap[info] = value;

      info  = "";
      value = "";
    }
  }

  std::ostream& operator << (std::ostream &out, const argInfoMap &m){
    std::map<std::string,std::string>::const_iterator it = m.iMap.begin();

    while(it != m.iMap.end()){
      out << it->first << " = " << it->second << '\n';
      ++it;
    }

    return out;
  }

  //---[ Kernel Info ]--------
  kernelInfo::kernelInfo() :
    occaKeywords(""),
    header(""),
    flags("") {}

  kernelInfo::kernelInfo(const kernelInfo &p) :
    occaKeywords(p.occaKeywords),
    header(p.header),
    flags(p.flags) {}

  kernelInfo& kernelInfo::operator = (const kernelInfo &p){
    occaKeywords = p.occaKeywords;
    header = p.header;
    flags  = p.flags;

    return *this;
  }

  kernelInfo& kernelInfo::operator += (const kernelInfo &p){
    header += p.header;
    flags  += p.flags;

    return *this;
  }

  std::string kernelInfo::salt() const {
    return (header + flags);
  }

  bool kernelInfo::isAnOccaDefine(const std::string &name){
    if((name == "OCCA_USING_CPU") ||
       (name == "OCCA_USING_GPU") ||

       (name == "OCCA_USING_SERIAL")   ||
       (name == "OCCA_USING_OPENMP")   ||
       (name == "OCCA_USING_OPENCL")   ||
       (name == "OCCA_USING_CUDA")     ||
       (name == "OCCA_USING_PTHREADS") ||
       (name == "OCCA_USING_COI")      ||

       (name == "occaInnerDim0") ||
       (name == "occaInnerDim1") ||
       (name == "occaInnerDim2") ||

       (name == "occaOuterDim0") ||
       (name == "occaOuterDim1") ||
       (name == "occaOuterDim2"))
      return true;

    return false;
  }

  void kernelInfo::addOCCAKeywords(const std::string &keywords){
    occaKeywords = keywords;
  }

  void kernelInfo::addIncludeDefine(const std::string &filename){
    header += "\n#include \"" + filename + "\"\n";
  }

  void kernelInfo::addInclude(const std::string &filename){
    header += '\n';
    header += readFile(filename);
    header += '\n';
  }

  void kernelInfo::removeDefine(const std::string &macro){
    if(!isAnOccaDefine(macro))
      header += "#undef " + macro + '\n';
  }

  void kernelInfo::addSource(const std::string &content){
    header += content;
  }

  void kernelInfo::addCompilerFlag(const std::string &f){
    flags += " " + f;
  }

  void kernelInfo::addCompilerIncludePath(const std::string &path){
#if (OCCA_OS & (LINUX_OS | OSX_OS))
    flags += " -I \"" + path + "\"";
#else
    flags += " /I \"" + path + "\"";
#endif
  }

  template <>
  void kernelInfo::addDefine(const std::string &macro, const std::string &value){
    std::stringstream ss;

    if(isAnOccaDefine(macro))
      ss << "#undef " << macro << "\n";

    // Make sure newlines are followed by escape characters
    std::string value2 = "";
    const int chars = value.size();

    for(int i = 0; i < chars; ++i){
      if(value[i] != '\n')
        value2 += value[i];
      else{
        if((i < (chars - 1))
           && (value[i] != '\\'))
          value2 += "\\\n";
        else
          value2 += '\n';
      }
    }

    if(value2[value2.size() - 1] != '\n')
      value2 += '\n';
    //==============

    ss << "#define " << macro << " " << value2 << '\n';

    header = ss.str() + header;
  }

  template <>
  void kernelInfo::addDefine(const std::string &macro, const float &value){
    std::stringstream ss;

    if(isAnOccaDefine(macro))
      ss << "#undef " << macro << "\n";

    ss << "#define " << macro << " ((float) " << std::setprecision(8) << value << ")\n";

    header = ss.str() + header;
  }

  template <>
  void kernelInfo::addDefine(const std::string &macro, const double &value){
    std::stringstream ss;

    if(isAnOccaDefine(macro))
      ss << "#undef " << macro << "\n";

    ss << "#define " << macro << " ((double) " << std::setprecision(16) << value << ")\n";

    header = ss.str() + header;
  }

  //---[ Device Info ]--------
  deviceInfo::deviceInfo(){}

  deviceInfo::deviceInfo(const deviceInfo &dInfo) :
    infos(dInfo.infos) {}

  deviceInfo& deviceInfo::operator = (const deviceInfo &dInfo){
    infos = dInfo.infos;

    return *this;
  }

  void deviceInfo::append(const std::string &key,
                          const std::string &value){
    if(infos.size() != 0)
      infos += ',';

    infos += key;
    infos += '=';
    infos += value;
  }

  //---[ Device :: Arg Info ]-
  argInfo::argInfo() :
    info(""),
    value("") {}

  argInfo::argInfo(const argInfo &ai) :
    info(ai.info),
    value(ai.value) {}

  argInfo& argInfo::operator = (const argInfo &ai){
    info  = ai.info;
    value = ai.value;

    return *this;
  }

  argInfo::argInfo(const std::string &info_) :
    info(info_),
    value("") {}

  argInfo::argInfo(const std::string &info_,
                   const std::string &value_) :
    info(info_),
    value(value_) {}

  const argInfo platformID("platformID");
  const argInfo deviceID("deviceID");

  const argInfo schedule("schedule");
  const argInfo chunk("chunk");

  const argInfo threadCount("threadCount");
  const argInfo pinnedCores("pinnedCores");
  //==========================
  //==================================

  //---[ Kernel ]---------------------
  kernel::kernel() :
    strMode(""),

    kHandle(NULL) {}
  kernel::kernel(kernel_v *kHandle_) :
    strMode( occa::modeToStr(kHandle_->mode()) ),

    kHandle(kHandle_) {}

  kernel::kernel(const kernel &k) :
    strMode(k.strMode),

    kHandle(k.kHandle) {}

  kernel& kernel::operator = (const kernel &k){
    strMode = k.strMode;

    kHandle = k.kHandle;

    return *this;
  }

  const std::string& kernel::mode(){
    return strMode;
  }

  kernel& kernel::buildFromSource(const std::string &filename,
                                  const std::string &functionName_,
                                  const kernelInfo &info_){
    kHandle->buildFromSource(filename, functionName_, info_);

    return *this;
  }

  kernel& kernel::buildFromBinary(const std::string &filename,
                                  const std::string &functionName_){
    kHandle->buildFromBinary(filename, functionName_);

    return *this;
  }

  kernel& kernel::loadFromLibrary(const char *cache,
                                  const std::string &functionName_){
    kHandle->loadFromLibrary(cache, functionName_);

    return *this;
  }

  void kernel::setWorkingDims(int dims, occa::dim inner, occa::dim outer){
    for(int i = 0; i < dims; ++i){
      inner[i] += (inner[i] ? 0 : 1);
      outer[i] += (outer[i] ? 0 : 1);
    }

    for(int i = dims; i < 3; ++i)
      inner[i] = outer[i] = 1;

    if(kHandle->nestedKernelCount == 0){
      kHandle->dims  = dims;
      kHandle->inner = inner;
      kHandle->outer = outer;
    }
    else{
      for(int k = 0; k < kHandle->nestedKernelCount; ++k)
        kHandle->nestedKernels[k].setWorkingDims(dims, inner, outer);
    }
  }

  int kernel::preferredDimSize(){
    OCCA_CHECK(kHandle->nestedKernelCount == 0,
               "Cannot get preferred size for fused kernels");

    return 1;
  }

  void kernel::clearArgumentList(){
    argumentCount = 0;
  }

  void kernel::addArgument(const int argPos,
                           const kernelArg &arg){
    if(argumentCount < (argPos + 1)){
      OCCA_CHECK(argPos < OCCA_MAX_ARGS,
                 "Kernels can only have at most [" << OCCA_MAX_ARGS << "] arguments,"
                 << " [" << argPos << "] arguments were set");

      argumentCount = (argPos + 1);
    }

    arguments[argPos] = arg;
  }

  void kernel::runFromArguments(){
    // [-] OCCA_MAX_ARGS = 25
#include "operators/occaRunFromArguments.cpp"

    return;
  }

#include "operators/occaOperatorDefinitions.cpp"

  void kernel::free(){
    if(kHandle->nestedKernelCount){
      for(int k = 0; k < kHandle->nestedKernelCount; ++k)
        kHandle->nestedKernels[k].free();

      delete [] kHandle->nestedKernels;
    }

    kHandle->free();

    delete kHandle;
  }

  kernelDatabase::kernelDatabase() :
    kernelName(""),
    modelKernelCount(0),
    kernelCount(0) {}

  kernelDatabase::kernelDatabase(const std::string kernelName_) :
    kernelName(kernelName_),
    modelKernelCount(0),
    kernelCount(0) {}

  kernelDatabase::kernelDatabase(const kernelDatabase &kdb) :
    kernelName(kdb.kernelName),

    modelKernelCount(kdb.modelKernelCount),
    modelKernelAvailable(kdb.modelKernelAvailable),

    kernelCount(kdb.kernelCount),
    kernels(kdb.kernels),
    kernelAllocated(kdb.kernelAllocated) {}


  kernelDatabase& kernelDatabase::operator = (const kernelDatabase &kdb){
    kernelName = kdb.kernelName;

    modelKernelCount     = kdb.modelKernelCount;
    modelKernelAvailable = kdb.modelKernelAvailable;

    kernelCount     = kdb.kernelCount;
    kernels         = kdb.kernels;
    kernelAllocated = kdb.kernelAllocated;

    return *this;
  }

  void kernelDatabase::modelKernelIsAvailable(const int id){
    OCCA_CHECK(0 <= id,
               "Model kernel for ID [" << id << "] was not found");

    if(modelKernelCount <= id){
      modelKernelCount = (id + 1);
      modelKernelAvailable.resize(modelKernelCount, false);
    }

    modelKernelAvailable[id] = true;
  }

  void kernelDatabase::addKernel(device d, kernel k){
    addKernel(d.dHandle->id_, k);
  }

  void kernelDatabase::addKernel(device_v *d, kernel k){
    addKernel(d->id_, k);
  }

  void kernelDatabase::addKernel(const int id, kernel k){
    OCCA_CHECK(0 <= id,
               "Model kernel for ID [" << id << "] was not found");

    if(kernelCount <= id){
      kernelCount = (id + 1);

      kernels.resize(kernelCount);
      kernelAllocated.resize(kernelCount, false);
    }

    kernels[id] = k;
    kernelAllocated[id] = true;
  }

  void kernelDatabase::loadKernelFromLibrary(device_v *d){
    addKernel(d, library::loadKernel(d, kernelName));
  }
  //==================================


  //---[ Memory ]---------------------
  memory::memory() :
    strMode(""),
    mHandle(NULL) {}

  memory::memory(void *uvaPtr){
    // Default to uvaPtr is actually a memory_v*
    memory_v *mHandle_ = (memory_v*) uvaPtr;

    ptrRangeMap_t::iterator it = uvaMap.find(uvaPtr);

    if(it != uvaMap.end())
      mHandle_ = it->second;

    strMode = occa::modeToStr(mHandle_->mode());
    mHandle = mHandle_;
  }

  memory::memory(memory_v *mHandle_) :
    strMode( occa::modeToStr(mHandle_->mode()) ),
    mHandle(mHandle_) {}

  memory::memory(const memory &m) :
    strMode(m.strMode),
    mHandle(m.mHandle) {}

  memory& memory::operator = (const memory &m){
    strMode = m.strMode;

    mHandle = m.mHandle;

    return *this;
  }

  const std::string& memory::mode(){
    return strMode;
  }

  void* memory::textureArg() const {
    return (void*) ((mHandle->textureInfo).arg);
  }

  device_v* memory::getOccaDeviceHandle(){
    return mHandle->dHandle;
  }

  memory_v* memory::getOccaMemoryHandle(){
    return mHandle;
  }

  void* memory::getMappedPointer(){
    return mHandle->mappedPtr;
  }

  void* memory::getMemoryHandle(){
    return mHandle->getMemoryHandle();
  }

  void* memory::getTextureHandle(){
    return mHandle->getTextureHandle();
  }

  void memory::placeInUva(){
    if( !(mHandle->dHandle->fakesUva()) ){
      mHandle->uvaPtr = mHandle->handle;
    }
    else if(mHandle->isMapped){
      mHandle->uvaPtr = mHandle->mappedPtr;
    }
    else{
      mHandle->uvaPtr = cpu::malloc(mHandle->size);
    }

    ptrRange_t uvaRange;

    uvaRange.start = (char*) (mHandle->uvaPtr);
    uvaRange.end   = (uvaRange.start + mHandle->size);

    uvaMap[uvaRange]                   = mHandle;
    mHandle->dHandle->uvaMap[uvaRange] = mHandle;

    // Needed for kernelArg.void_ -> mHandle checks
    if(mHandle->uvaPtr != mHandle->handle)
      uvaMap[mHandle->handle] = mHandle;
  }

  void memory::manage(){
    placeInUva();

    mHandle->isManaged = true;
  }

  void memory::syncToDevice(const uintptr_t bytes){
    if(mHandle->dHandle->fakesUva()){
      uintptr_t bytes_ = ((bytes == 0) ? mHandle->size : bytes);

      copyTo(mHandle->uvaPtr, bytes_);

      mHandle->uva_inDevice = true;
      mHandle->uva_isDirty  = false;

      removeFromDirtyMap(mHandle);
    }
  }

  void memory::syncFromDevice(const uintptr_t bytes){
    if(mHandle->dHandle->fakesUva()){
      uintptr_t bytes_ = ((bytes == 0) ? mHandle->size : bytes);

      copyFrom(mHandle->uvaPtr, bytes_);

      mHandle->uva_inDevice = false;
      mHandle->uva_isDirty  = false;

      removeFromDirtyMap(mHandle);
    }
  }

  bool memory::uvaIsDirty(){
    return (mHandle && (mHandle->uva_isDirty));
  }

  void memory::uvaMarkDirty(){
    if(mHandle)
      mHandle->uva_isDirty = true;
  }

  void memory::uvaMarkClean(){
    if(mHandle)
      mHandle->uva_isDirty = false;
  }

  void memory::copyFrom(const void *src,
                        const uintptr_t bytes,
                        const uintptr_t offset){

    mHandle->copyFrom(src, bytes, offset);
  }

  void memory::copyFrom(const memory &src,
                        const uintptr_t bytes,
                        const uintptr_t destOffset,
                        const uintptr_t srcOffset){

    if(mHandle->dHandle == src.mHandle->dHandle){
      mHandle->copyFrom(src.mHandle, bytes, destOffset, srcOffset);
    }
    else{
      memory_v *srcHandle  = src.mHandle;
      memory_v *destHandle = mHandle;

      const occa::mode modeS = srcHandle->mode();
      const occa::mode modeD = destHandle->mode();

      if(modeS & onChipModes){
        destHandle->copyFrom(srcHandle->getMemoryHandle(),
                             bytes, destOffset);
      }
      else if(modeD & onChipModes){
        srcHandle->copyTo(destHandle->getMemoryHandle(),
                          bytes, srcOffset);
      }
      else{
        OCCA_CHECK(((modeS == CUDA) && (modeD == CUDA)),
                   "Peer-to-peer is not supported between ["
                   << modeToStr(modeS) << "] and ["
                   << modeToStr(modeD) << "]");

#if OCCA_CUDA_ENABLED
        CUDADeviceData_t &srcDevData  =
          *((CUDADeviceData_t*) srcHandle->dHandle->data);

        CUDADeviceData_t &destDevData =
          *((CUDADeviceData_t*) destHandle->dHandle->data);

        CUdeviceptr srcMem  = *(((CUdeviceptr*) srcHandle->handle)  + srcOffset);
        CUdeviceptr destMem = *(((CUdeviceptr*) destHandle->handle) + destOffset);

        if(!srcDevData.p2pEnabled)
          cuda::enablePeerToPeer(srcDevData.context);

        if(!destDevData.p2pEnabled)
          cuda::enablePeerToPeer(destDevData.context);

        cuda::checkPeerToPeer(destDevData.device,
                              srcDevData.device);

        cuda::peerToPeerMemcpy(destDevData.device,
                               destDevData.context,
                               destMem,

                               srcDevData.device,
                               srcDevData.context,
                               srcMem,

                               bytes,
                               (CUstream) srcHandle->dHandle->currentStream);
#endif
      }
    }
  }

  void memory::copyTo(void *dest,
                      const uintptr_t bytes,
                      const uintptr_t offset){

    mHandle->copyTo(dest, bytes, offset);
  }

  void memory::copyTo(memory &dest,
                      const uintptr_t bytes,
                      const uintptr_t destOffset,
                      const uintptr_t srcOffset){

    if(mHandle->dHandle == dest.mHandle->dHandle){
      mHandle->copyTo(dest.mHandle, bytes, destOffset, srcOffset);
    }
    else{
      memory_v *srcHandle  = mHandle;
      memory_v *destHandle = dest.mHandle;

      const occa::mode modeS = srcHandle->mode();
      const occa::mode modeD = destHandle->mode();

      if(modeS & onChipModes){
        destHandle->copyFrom(srcHandle->getMemoryHandle(),
                             bytes, srcOffset);
      }
      else if(modeD & onChipModes){
        srcHandle->copyTo(destHandle->getMemoryHandle(),
                          bytes, destOffset);
      }
      else{
        OCCA_CHECK(((modeS == CUDA) && (modeD == CUDA)),
                   "Peer-to-peer is not supported between ["
                   << modeToStr(modeS) << "] and ["
                   << modeToStr(modeD) << "]");

#if OCCA_CUDA_ENABLED
        CUDADeviceData_t &srcDevData  =
          *((CUDADeviceData_t*) srcHandle->dHandle->data);

        CUDADeviceData_t &destDevData =
          *((CUDADeviceData_t*) destHandle->dHandle->data);

        CUdeviceptr srcMem  = *(((CUdeviceptr*) srcHandle->handle)  + srcOffset);
        CUdeviceptr destMem = *(((CUdeviceptr*) destHandle->handle) + destOffset);

        cuda::peerToPeerMemcpy(destDevData.device,
                               destDevData.context,
                               destMem,

                               srcDevData.device,
                               srcDevData.context,
                               srcMem,

                               bytes,
                               (CUstream) srcHandle->dHandle->currentStream);
#endif
      }
    }
  }

  void memory::asyncCopyFrom(const void *src,
                             const uintptr_t bytes,
                             const uintptr_t offset){

    mHandle->asyncCopyFrom(src, bytes, offset);
  }

  void memory::asyncCopyFrom(const memory &src,
                             const uintptr_t bytes,
                             const uintptr_t destOffset,
                             const uintptr_t srcOffset){

    if(mHandle->dHandle == src.mHandle->dHandle){
      mHandle->asyncCopyFrom(src.mHandle, bytes, destOffset, srcOffset);
    }
    else{
      memory_v *srcHandle  = src.mHandle;
      memory_v *destHandle = mHandle;

      const occa::mode modeS = srcHandle->mode();
      const occa::mode modeD = destHandle->mode();

      if(modeS & onChipModes){
        destHandle->asyncCopyFrom(srcHandle->getMemoryHandle(),
                             bytes, destOffset);
      }
      else if(modeD & onChipModes){
        srcHandle->asyncCopyTo(destHandle->getMemoryHandle(),
                          bytes, srcOffset);
      }
      else{
        OCCA_CHECK(((modeS == CUDA) && (modeD == CUDA)),
                   "Peer-to-peer is not supported between ["
                   << modeToStr(modeS) << "] and ["
                   << modeToStr(modeD) << "]");

#if OCCA_CUDA_ENABLED
        CUDADeviceData_t &srcDevData  =
          *((CUDADeviceData_t*) srcHandle->dHandle->data);

        CUDADeviceData_t &destDevData =
          *((CUDADeviceData_t*) destHandle->dHandle->data);

        CUdeviceptr srcMem  = *(((CUdeviceptr*) srcHandle->handle)  + srcOffset);
        CUdeviceptr destMem = *(((CUdeviceptr*) destHandle->handle) + destOffset);

        cuda::asyncPeerToPeerMemcpy(destDevData.device,
                                    destDevData.context,
                                    destMem,

                                    srcDevData.device,
                                    srcDevData.context,
                                    srcMem,

                                    bytes,
                                    (CUstream) srcHandle->dHandle->currentStream);
#endif
      }
    }
  }

  void memory::asyncCopyTo(void *dest,
                           const uintptr_t bytes,
                           const uintptr_t offset){

    mHandle->asyncCopyTo(dest, bytes, offset);
  }

  void memory::asyncCopyTo(memory &dest,
                           const uintptr_t bytes,
                           const uintptr_t destOffset,
                           const uintptr_t srcOffset){

    if(mHandle->dHandle == dest.mHandle->dHandle){
      mHandle->asyncCopyTo(dest.mHandle, bytes, destOffset, srcOffset);
    }
    else{
      memory_v *srcHandle  = mHandle;
      memory_v *destHandle = dest.mHandle;

      const occa::mode modeS = srcHandle->mode();
      const occa::mode modeD = destHandle->mode();

      if(modeS & onChipModes){
        destHandle->asyncCopyFrom(srcHandle->getMemoryHandle(),
                                  bytes, destOffset);
      }
      else if(modeD & onChipModes){
        srcHandle->asyncCopyTo(destHandle->getMemoryHandle(),
                               bytes, srcOffset);
      }
      else{
        OCCA_CHECK(((modeS == CUDA) && (modeD == CUDA)),
                   "Peer-to-peer is not supported between ["
                   << modeToStr(modeS) << "] and ["
                   << modeToStr(modeD) << "]");

#if OCCA_CUDA_ENABLED
        CUDADeviceData_t &srcDevData  =
          *((CUDADeviceData_t*) srcHandle->dHandle->data);

        CUDADeviceData_t &destDevData =
          *((CUDADeviceData_t*) destHandle->dHandle->data);

        CUdeviceptr srcMem  = *(((CUdeviceptr*) srcHandle->handle)  + srcOffset);
        CUdeviceptr destMem = *(((CUdeviceptr*) destHandle->handle) + destOffset);

        cuda::asyncPeerToPeerMemcpy(destDevData.device,
                                    destDevData.context,
                                    destMem,

                                    srcDevData.device,
                                    srcDevData.context,
                                    srcMem,

                                    bytes,
                                    (CUstream) srcHandle->dHandle->currentStream);
#endif
      }
    }
  }

  void memcpy(void *dest, void *src,
              const uintptr_t bytes,
              const int flags){

    memcpy(dest, src, bytes, flags, false);
  }

  void asyncMemcpy(void *dest, void *src,
                   const uintptr_t bytes,
                   const int flags){

    memcpy(dest, src, bytes, flags, true);
  }

  void memcpy(void *dest, void *src,
              const uintptr_t bytes,
              const int flags,
              const bool isAsync){

    ptrRangeMap_t::iterator srcIt  = uvaMap.end();
    ptrRangeMap_t::iterator destIt = uvaMap.end();

    if(flags & occa::autoDetect){
      srcIt  = uvaMap.find(src);
      destIt = uvaMap.find(dest);
    }
    else{
      if(flags & srcInUva)
        srcIt  = uvaMap.find(src);

      if(flags & destInUva)
        destIt  = uvaMap.find(dest);
    }

    occa::memory_v *srcMem  = ((srcIt != uvaMap.end())  ? (srcIt->second)  : NULL);
    occa::memory_v *destMem = ((destIt != uvaMap.end()) ? (destIt->second) : NULL);

    const uintptr_t srcOff  = (srcMem  ? (((char*) src)  - ((char*) srcMem->uvaPtr))  : 0);
    const uintptr_t destOff = (destMem ? (((char*) dest) - ((char*) destMem->uvaPtr)) : 0);

    const bool usingSrcPtr  = ((srcMem  == NULL) || (srcMem->isManaged));
    const bool usingDestPtr = ((destMem == NULL) || (destMem->isManaged));

    if(usingSrcPtr && usingDestPtr){
      ::memcpy(dest, src, bytes);
    }
    else if(usingSrcPtr){
      if(!isAsync)
        destMem->copyFrom(src, bytes, destOff);
      else
        destMem->asyncCopyFrom(src, bytes, destOff);
    }
    else if(usingDestPtr){
      if(!isAsync)
        srcMem->copyTo(dest, bytes, srcOff);
      else
        srcMem->asyncCopyTo(dest, bytes, srcOff);
    }
    else {
      // Auto-detects peer-to-peer stuff
      occa::memory srcMemory(srcMem);
      occa::memory destMemory(destMem);

      if(!isAsync)
        srcMemory.copyTo(destMemory, bytes, destOff, srcOff);
      else
        srcMemory.asyncCopyTo(destMemory, bytes, destOff, srcOff);
    }
  }

  void memcpy(memory &dest,
              const void *src,
              const uintptr_t bytes,
              const uintptr_t offset){

    dest.copyFrom(src, bytes, offset);
  }

  void memcpy(void *dest,
              memory &src,
              const uintptr_t bytes,
              const uintptr_t offset){

    src.copyTo(dest, bytes, offset);
  }

  void memcpy(memory &dest,
              memory &src,
              const uintptr_t bytes,
              const uintptr_t destOffset,
              const uintptr_t srcOffset){

    src.copyTo(dest, bytes, destOffset, srcOffset);
  }

  void asyncMemcpy(memory &dest,
                   const void *src,
                   const uintptr_t bytes,
                   const uintptr_t offset){

    dest.asyncCopyFrom(src, bytes, offset);
  }

  void asyncMemcpy(void *dest,
                   memory &src,
                   const uintptr_t bytes,
                   const uintptr_t offset){

    src.asyncCopyTo(dest, bytes, offset);
  }

  void asyncMemcpy(memory &dest,
                   memory &src,
                   const uintptr_t bytes,
                   const uintptr_t destOffset,
                   const uintptr_t srcOffset){

    src.asyncCopyTo(dest, bytes, destOffset, srcOffset);
  }

  void memory::swap(memory &m){
    std::string strMode2 = m.strMode;
    m.strMode            = strMode;
    strMode              = strMode2;

    memory_v *mHandle2 = m.mHandle;
    m.mHandle          = mHandle;
    mHandle            = mHandle2;
  }

  void memory::free(){
    mHandle->dHandle->bytesAllocated -= (mHandle->size);

    if(mHandle->uvaPtr){
      uvaMap.erase(mHandle->uvaPtr);
      mHandle->dHandle->uvaMap.erase(mHandle->uvaPtr);

      // CPU case where memory is shared
      if(mHandle->uvaPtr != mHandle->handle){
        uvaMap.erase(mHandle->handle);
        mHandle->dHandle->uvaMap.erase(mHandle->uvaPtr);

        ::free(mHandle->uvaPtr);
      }
    }

    if( !(mHandle->isMapped) )
      mHandle->free();
    else
      mHandle->mappedFree();

    delete mHandle;
  }
  //==================================


  //---[ Device ]---------------------
  device::device() :
    dHandle(NULL) {}

  device::device(device_v *dHandle_) :
    dHandle(dHandle_) {}

  device::device(deviceInfo &dInfo){
    setup(dInfo);
  }

  device::device(const std::string &infos){
    setup(infos);
  }

  device::device(const device &d) :
    dHandle(d.dHandle) {}

  device& device::operator = (const device &d){
    dHandle = d.dHandle;

    return *this;
  }

  void device::setupHandle(occa::mode m){
    switch(m){

    case Serial:{
      dHandle = new device_t<Serial>();
      break;
    }
    case OpenMP:{
#if OCCA_OPENMP_ENABLED
      dHandle = new device_t<OpenMP>();
#else
      std::cout << "OCCA mode [OpenMP] is not enabled, defaulting to [Serial] mode\n";
      dHandle = new device_t<Serial>();
#endif
      break;
    }
    case OpenCL:{
#if OCCA_OPENCL_ENABLED
      dHandle = new device_t<OpenCL>();
#else
      std::cout << "OCCA mode [OpenCL] is not enabled, defaulting to [Serial] mode\n";
      dHandle = new device_t<Serial>();
#endif
      break;
    }
    case CUDA:{
#if OCCA_CUDA_ENABLED
      dHandle = new device_t<CUDA>();
#else
      std::cout << "OCCA mode [CUDA] is not enabled, defaulting to [Serial] mode\n";
      dHandle = new device_t<Serial>();
#endif
      break;
    }
    case Pthreads:{
      std::cout << "OCCA mode [Pthreads] is still in development-mode (unstable)\n";
      dHandle = new device_t<Pthreads>();
      break;
    }
    case COI:{
#if OCCA_COI_ENABLED
      std::cout << "OCCA mode [COI] is deprecated (unstable)\n";
      dHandle = new device_t<COI>();
#else
      std::cout << "OCCA mode [COI] is not enabled, defaulting to [Serial] mode\n";
      dHandle = new device_t<Serial>();
#endif
      break;
    }
    default:{
      std::cout << "Unsupported OCCA mode given, defaulting to [Serial] mode\n";
      dHandle = new device_t<Serial>();
    }
    }
  }

  void device::setupHandle(const std::string &m){
    setupHandle( strToMode(m) );
  }

  void device::setup(deviceInfo &dInfo){
    setup(dInfo.infos);
  }

  void device::setup(const std::string &infos){
    argInfoMap aim(infos);

    OCCA_CHECK(aim.has("mode"),
               "OCCA mode not given");

    // Load [mode] from aim
    occa::mode m = strToMode(aim.get("mode"));

    setupHandle(m);

    dHandle->setup(aim);

    dHandle->modelID_ = library::deviceModelID(dHandle->getIdentifier());
    dHandle->id_      = library::genDeviceID();

    if(aim.has("UVA")){
      if(upStringCheck(aim.get("UVA"), "enabled"))
        dHandle->uvaEnabled_ = true;
      else
        dHandle->uvaEnabled_ = false;
    }
    else
      dHandle->uvaEnabled_ = uvaEnabledByDefault_f;

    dHandle->currentStream = createStream();
  }

  void device::setup(occa::mode m,
                     const int arg1, const int arg2){
    setupHandle(m);

    argInfoMap aim;

    switch(m){
    case Serial:{
      // Do Nothing
      break;
    }
    case OpenMP:{
      // Do Nothing, maybe add thread order next, dynamic static, etc
      break;
    }
    case OpenCL:{
      aim.set("platformID", arg1);
      aim.set("deviceID"  , arg2);
      break;
    }
    case CUDA:{
      aim.set("deviceID", arg1);
      break;
    }
    case Pthreads:{
      aim.set("threadCount", arg1);
      aim.set("pinningInfo", arg2);
      break;
    }
    case COI:{
      aim.set("deviceID", arg1);
      break;
    }
    }

    dHandle->setup(aim);

    dHandle->modelID_ = library::deviceModelID(dHandle->getIdentifier());
    dHandle->id_      = library::genDeviceID();

    dHandle->currentStream = createStream();
  }

  void device::setup(occa::mode m,
                     const argInfo &arg1){
    setupHandle(m);

    argInfoMap aim;

    aim.set(arg1.info, arg1.value);

    dHandle->setup(aim);

    dHandle->modelID_ = library::deviceModelID(dHandle->getIdentifier());
    dHandle->id_      = library::genDeviceID();

    dHandle->currentStream = createStream();
  }

  void device::setup(occa::mode m,
                     const argInfo &arg1, const argInfo &arg2){
    setupHandle(m);

    argInfoMap aim;

    aim.set(arg1.info, arg1.value);
    aim.set(arg2.info, arg2.value);

    dHandle->setup(aim);

    dHandle->modelID_ = library::deviceModelID(dHandle->getIdentifier());
    dHandle->id_      = library::genDeviceID();

    dHandle->currentStream = createStream();
  }

  void device::setup(occa::mode m,
                     const argInfo &arg1, const argInfo &arg2, const argInfo &arg3){
    setupHandle(m);

    argInfoMap aim;

    aim.set(arg1.info, arg1.value);
    aim.set(arg2.info, arg2.value);
    aim.set(arg3.info, arg3.value);

    dHandle->setup(aim);

    dHandle->modelID_ = library::deviceModelID(dHandle->getIdentifier());
    dHandle->id_      = library::genDeviceID();

    dHandle->currentStream = createStream();
  }


  void device::setup(const std::string &m,
                     const int arg1, const int arg2){
    setup(strToMode(m), arg1, arg2);
  }

  void device::setup(const std::string &m,
                     const argInfo &arg1){
    setup(strToMode(m), arg1);
  }

  void device::setup(const std::string &m,
                     const argInfo &arg1, const argInfo &arg2){
    setup(strToMode(m), arg1, arg2);
  }

  void device::setup(const std::string &m,
                     const argInfo &arg1, const argInfo &arg2, const argInfo &arg3){
    setup(strToMode(m), arg1, arg2, arg3);
  }

  uintptr_t device::bytesAllocated() const {
    return dHandle->bytesAllocated;
  }

  deviceIdentifier device::getIdentifier() const {
    return dHandle->getIdentifier();
  }

  void device::setCompiler(const std::string &compiler_){
    dHandle->setCompiler(compiler_);
  }

  void device::setCompilerEnvScript(const std::string &compilerEnvScript_){
    dHandle->setCompilerEnvScript(compilerEnvScript_);
  }

  void device::setCompilerFlags(const std::string &compilerFlags_){
    dHandle->setCompilerFlags(compilerFlags_);
  }

  std::string& device::getCompiler(){
    return dHandle->compiler;
  }

  std::string& device::getCompilerEnvScript(){
    return dHandle->compilerEnvScript;
  }

  std::string& device::getCompilerFlags(){
    return dHandle->compilerFlags;
  }

  int device::modelID(){
    return dHandle->modelID_;
  }

  int device::id(){
    return dHandle->id_;
  }

  int device::modeID(){
    return dHandle->mode();
  }

  const std::string& device::mode(){
    return strMode;
  }

  void device::flush(){
    dHandle->flush();
  }

  void device::finish(){
    if(dHandle->fakesUva()){
      const size_t dirtyEntries = uvaDirtyMemory.size();

      if(dirtyEntries){
        for(size_t i = 0; i < dirtyEntries; ++i){
          occa::memory_v *mem = uvaDirtyMemory[i];

          mem->asyncCopyTo(mem->uvaPtr);

          mem->uva_inDevice = false;
          mem->uva_isDirty  = false;
        }

        uvaDirtyMemory.clear();
      }
    }

    dHandle->finish();
  }

  void device::waitFor(streamTag tag){
    dHandle->waitFor(tag);
  }

  stream device::createStream(){
    dHandle->streams.push_back( dHandle->createStream() );
    return dHandle->streams.back();
  }

  stream device::getStream(){
    return dHandle->currentStream;
  }

  void device::setStream(stream s){
    dHandle->currentStream = s;
  }

  stream device::wrapStream(void *handle_){
    return dHandle->wrapStream(handle_);
  }

  streamTag device::tagStream(){
    return dHandle->tagStream();
  }

  double device::timeBetween(const streamTag &startTag, const streamTag &endTag){
    return dHandle->timeBetween(startTag, endTag);
  }

  void device::freeStream(stream s){
    const int streamCount = dHandle->streams.size();

    for(int i = 0; i < streamCount; ++i){
      if(dHandle->streams[i] == s){
        dHandle->freeStream(dHandle->streams[i]);
        dHandle->streams.erase(dHandle->streams.begin() + i);

        break;
      }
    }
  }

  kernel device::buildKernel(const std::string &str,
                             const std::string &functionName,
                             const kernelInfo &info_){

    if(fileExists(str))
      return buildKernelFromSource(str, functionName, info_);
    else
      return buildKernelFromString(str, functionName, info_);
  }

  kernel device::buildKernelFromString(const std::string &content,
                                       const std::string &functionName,
                                       const int language){

    return buildKernelFromString(content, functionName, defaultKernelInfo, language);
  }

  kernel device::buildKernelFromString(const std::string &content,
                                       const std::string &functionName,
                                       const kernelInfo &info_,
                                       const int language){

    kernelInfo info = info_;

    dHandle->addOccaHeadersToInfo(info);

    const std::string cachedBinary = getContentCachedName(content, dHandle->getInfoSalt(info));

    std::string prefix, cacheName;

    getFilePrefixAndName(cachedBinary, prefix, cacheName);

    std::string h_cacheName = prefix + "h_" + cacheName;

    if(language & occa::usingOKL)
      h_cacheName += ".okl";
    else if(language & occa::usingOFL)
      h_cacheName += ".ofl";
    else
      h_cacheName += ".occa";

    if(!haveFile(h_cacheName)){
      waitForFile(h_cacheName);
      waitForFile(cachedBinary);

      return buildKernelFromBinary(cachedBinary, functionName);
    }

    writeToFile(h_cacheName, content);
    releaseFile(h_cacheName);

    return buildKernelFromSource(h_cacheName, functionName, info_);
  }

  kernel device::buildKernelFromSource(const std::string &filename,
                                       const std::string &functionName,
                                       const kernelInfo &info_){

    const bool usingParser = fileNeedsParser(filename);

    kernel ker;

    ker.strMode = strMode;

    kernel_v *&k = ker.kHandle;

    if(usingParser){
      k          = new kernel_t<Serial>;
      k->dHandle = new device_t<Serial>();

      kernelInfo info = info_;

      const std::string cachedBinary  = k->getCachedBinaryName(filename, info);
      const std::string iCachedBinary = getMidCachedBinaryName(cachedBinary, "i");

      k->metaInfo = parseFileForFunction(filename,
                                         cachedBinary,
                                         functionName,
                                         info_);

      info = defaultKernelInfo;
      info.addDefine("OCCA_LAUNCH_KERNEL", 1);

      struct stat buffer;
      bool fileExists = (stat(cachedBinary.c_str(), &buffer) == 0);

      if(fileExists){
        if(verboseCompilation_f)
          std::cout << "Found cached binary of [" << filename << "] in [" << cachedBinary << "]\n";

        k->buildFromBinary(cachedBinary, functionName);
      }
      else
        k->buildFromSource(iCachedBinary, functionName, info);

      k->nestedKernelCount = k->metaInfo.nestedKernels;

      std::stringstream ss;
      k->nestedKernels = new kernel[k->metaInfo.nestedKernels];

      for(int ki = 0; ki < k->metaInfo.nestedKernels; ++ki){
        ss << ki;

        const std::string sKerName = k->metaInfo.baseName + ss.str();

        ss.str("");

        kernel &sKer = k->nestedKernels[ki];

        sKer.strMode = strMode;

        sKer.kHandle = dHandle->buildKernelFromSource(iCachedBinary,
                                                      sKerName,
                                                      info_);

        sKer.kHandle->metaInfo               = k->metaInfo;
        sKer.kHandle->metaInfo.name          = sKerName;
        sKer.kHandle->metaInfo.nestedKernels = 0;
        sKer.kHandle->metaInfo.removeArg(0); // remove nestedKernels **
      }
    }
    else{
      k          = dHandle->buildKernelFromSource(filename, functionName, info_);
      k->dHandle = dHandle;
    }

    return ker;
  }

  kernel device::buildKernelFromBinary(const std::string &filename,
                                       const std::string &functionName){
    kernel ker;

    ker.strMode = strMode;

    ker.kHandle          = dHandle->buildKernelFromBinary(filename, functionName);
    ker.kHandle->dHandle = dHandle;

    return ker;
  }

  void device::cacheKernelInLibrary(const std::string &filename,
                                    const std::string &functionName,
                                    const kernelInfo &info_){
    dHandle->cacheKernelInLibrary(filename, functionName, info_);
  }

  kernel device::loadKernelFromLibrary(const char *cache,
                                       const std::string &functionName){
    kernel ker;

    ker.strMode = strMode;

    ker.kHandle          = dHandle->loadKernelFromLibrary(cache, functionName);
    ker.kHandle->dHandle = dHandle;

    return ker;
  }

  kernel device::buildKernelFromLoopy(const std::string &filename,
                                      const std::string &functionName,
                                      const int useLoopyOrFloopy){

    return buildKernelFromLoopy(filename,
                                functionName,
                                defaultKernelInfo,
                                useLoopyOrFloopy);
  }

  kernel device::buildKernelFromLoopy(const std::string &filename,
                                      const std::string &functionName,
                                      const kernelInfo &info_,
                                      const int useLoopyOrFloopy){

    std::string cachedBinary = getCachedName(filename, dHandle->getInfoSalt(info_));

    struct stat buffer;
    bool fileExists = (stat(cachedBinary.c_str(), &buffer) == 0);

    if(fileExists){
      if(verboseCompilation_f)
        std::cout << "Found loo.py cached binary of [" << filename << "] in [" << cachedBinary << "]\n";

      return buildKernelFromBinary(cachedBinary, functionName);
    }

    std::string prefix, cacheName;

    getFilePrefixAndName(cachedBinary, prefix, cacheName);

    const std::string defsFile = prefix + "loopy1_" + cacheName + ".defs";
    const std::string clFile = prefix + "loopy2_" + cacheName + ".ocl";

    writeToFile(defsFile, info_.header);

    std::string loopyLang = "loopy";

    if(useLoopyOrFloopy == occa::useFloopy)
      loopyLang = "fpp";

    std::stringstream command;

    command << "loopy --lang=" << loopyLang
            << " --occa-defines=" << defsFile << ' '
            << " --occa-add-dummy-arg "
            << filename << ' ' << clFile;

    const std::string &sCommand = command.str();

    if(verboseCompilation_f)
      std::cout << sCommand << '\n';

#if (OCCA_OS & (LINUX_OS | OSX_OS))
    const int compileError = system(sCommand.c_str());
#else
    const int compileError = system(("\"" +  sCommand + "\"").c_str());
#endif

    if(compileError)
      OCCA_CHECK(false, "Compilation error");

    return buildKernelFromSource(clFile, functionName);
  }

  memory device::wrapMemory(void *handle_,
                            const uintptr_t bytes){
    memory mem;

    mem.strMode = strMode;

    mem.mHandle = dHandle->wrapMemory(handle_, bytes);
    mem.mHandle->dHandle = dHandle;

    return mem;
  }

  void* device::wrapManagedMemory(void *handle_,
                                  const uintptr_t bytes){
    memory mem = wrapMemory(handle_, bytes);

    mem.manage();

    return mem.mHandle->uvaPtr;
  }

  memory device::wrapTexture(void *handle_,
                             const int dim, const occa::dim &dims,
                             occa::formatType type, const int permissions){

    OCCA_CHECK((dim == 1) || (dim == 2),
               "Textures of [" << dim << "D] are not supported,"
               << "only 1D or 2D are supported at the moment");

    memory mem;

    mem.strMode = strMode;

    mem.mHandle = dHandle->wrapTexture(handle_,
                                       dim, dims,
                                       type, permissions);
    mem.mHandle->dHandle = dHandle;

    return mem;
  }

  void* device::wrapManagedTexture(void *handle_,
                                   const int dim, const occa::dim &dims,
                                   occa::formatType type, const int permissions){

    memory mem = wrapTexture(handle_, dim, dims, type, permissions);

    mem.manage();

    return mem.mHandle->uvaPtr;
  }

  memory device::malloc(const uintptr_t bytes,
                        void *src){
    memory mem;

    mem.strMode = strMode;

    mem.mHandle          = dHandle->malloc(bytes, src);
    mem.mHandle->dHandle = dHandle;

    dHandle->bytesAllocated += bytes;

    return mem;
  }

  void* device::managedAlloc(const uintptr_t bytes,
                             void *src){
    memory mem = malloc(bytes, src);

    mem.manage();

    return mem.mHandle->uvaPtr;
  }

  void* device::uvaAlloc(const uintptr_t bytes,
                         void *src){
    memory mem = malloc(bytes, src);

    mem.placeInUva();

    return mem.mHandle->uvaPtr;
  }

  void* device::managedUvaAlloc(const uintptr_t bytes,
                                void *src){
    memory mem = malloc(bytes, src);

    mem.manage();

    return mem.mHandle->uvaPtr;
  }

  memory device::textureAlloc(const int dim, const occa::dim &dims,
                              void *src,
                              occa::formatType type, const int permissions){
    OCCA_CHECK((dim == 1) || (dim == 2),
               "Textures of [" << dim << "D] are not supported,"
               << "only 1D or 2D are supported at the moment");

    OCCA_CHECK(src != NULL,
               "Non-NULL source is required for [textureAlloc] (texture allocation)");

    memory mem;

    mem.strMode = strMode;

    mem.mHandle      = dHandle->textureAlloc(dim, dims, src, type, permissions);
    mem.mHandle->dHandle = dHandle;

    dHandle->bytesAllocated += (type.bytes() *
                                ((dim == 2) ?
                                 (dims[0] * dims[1]) :
                                 (dims[0]          )));

    return mem;
  }

  void* device::managedTextureAlloc(const int dim, const occa::dim &dims,
                                    void *src,
                                    occa::formatType type, const int permissions){

    memory mem = textureAlloc(dim, dims, src, type, permissions);

    mem.manage();

    return mem.mHandle->uvaPtr;
  }

  memory device::mappedAlloc(const uintptr_t bytes,
                             void *src){
    memory mem;

    mem.strMode = strMode;

    mem.mHandle          = dHandle->mappedAlloc(bytes, src);
    mem.mHandle->dHandle = dHandle;

    dHandle->bytesAllocated += bytes;

    return mem;
  }

  void* device::managedMappedAlloc(const uintptr_t bytes,
                                   void *src){
    memory mem = mappedAlloc(bytes, src);

    mem.manage();

    return mem.mHandle->uvaPtr;
  }

  void device::free(){
    const int streamCount = dHandle->streams.size();

    for(int i = 0; i < streamCount; ++i)
      dHandle->freeStream(dHandle->streams[i]);

    dHandle->free();

    delete dHandle;
  }

  int device::simdWidth(){
    return dHandle->simdWidth();
  }

  mutex_t deviceListMutex;
  std::vector<device> deviceList;

  std::vector<device>& getDeviceList(){

    deviceListMutex.lock();

    if(deviceList.size()){
      deviceListMutex.unlock();
      return deviceList;
    }

    device_t<Serial>::appendAvailableDevices(deviceList);

#if OCCA_OPENMP_ENABLED
    device_t<OpenMP>::appendAvailableDevices(deviceList);
#endif
#if OCCA_PTHREADS_ENABLED
    device_t<Pthreads>::appendAvailableDevices(deviceList);
#endif
#if OCCA_OPENCL_ENABLED
    device_t<OpenCL>::appendAvailableDevices(deviceList);
#endif
#if OCCA_CUDA_ENABLED
    device_t<CUDA>::appendAvailableDevices(deviceList);
#endif
#if OCCA_COI_ENABLED
    device_t<COI>::appendAvailableDevices(deviceList);
#endif

    deviceListMutex.unlock();

    return deviceList;
  }
  void printAvailableDevices(){
    std::stringstream ss;
    ss << "==============o=======================o==========================================\n";
    ss << cpu::getDeviceListInfo();
#if OCCA_OPENCL_ENABLED
    ss << "==============o=======================o==========================================\n";
    ss << cl::getDeviceListInfo();
#endif
#if OCCA_CUDA_ENABLED
    ss << "==============o=======================o==========================================\n";
    ss << cuda::getDeviceListInfo();
#endif
#if OCCA_COI_ENABLED
    ss << "==============o=======================o==========================================\n";
    ss << coi::getDeviceListInfo();
#endif
    ss << "==============o=======================o==========================================\n";

    std::cout << ss.str();
  }

  deviceIdentifier::deviceIdentifier() :
    mode_(Serial) {}

  deviceIdentifier::deviceIdentifier(occa::mode m,
                                     const char *c, const size_t chars){
    mode_ = m;
    load(c, chars);
  }

  deviceIdentifier::deviceIdentifier(occa::mode m, const std::string &s){
    mode_ = m;
    load(s);
  }

  deviceIdentifier::deviceIdentifier(const deviceIdentifier &di) :
    mode_(di.mode_),
    flagMap(di.flagMap) {}

  deviceIdentifier& deviceIdentifier::operator = (const deviceIdentifier &di){
    mode_ = di.mode_;
    flagMap = di.flagMap;

    return *this;
  }

  void deviceIdentifier::load(const char *c, const size_t chars){
    const char *c1 = c;

    while((c1 < (c + chars)) && (*c1 != '\0')){
      const char *c2 = c1;
      const char *c3;

      while(*c2 != '|')
        ++c2;

      c3 = (c2 + 1);

      while((c3 < (c + chars)) &&
            (*c3 != '\0') && (*c3 != '|'))
        ++c3;

      flagMap[std::string(c1, c2 - c1)] = std::string(c2 + 1, c3 - c2 - 1);

      c1 = (c3 + 1);
    }
  }

  void deviceIdentifier::load(const std::string &s){
    return load(s.c_str(), s.size());
  }

  std::string deviceIdentifier::flattenFlagMap() const {
    std::string ret = "";

    cFlagMapIterator it = flagMap.begin();

    if(it == flagMap.end())
      return "";

    ret += it->first;
    ret += '|';
    ret += it->second;
    ++it;

    while(it != flagMap.end()){
      ret += '|';
      ret += it->first;
      ret += '|';
      ret += it->second;

      ++it;
    }

    return ret;
  }

  int deviceIdentifier::compare(const deviceIdentifier &b) const {
    if(mode_ != b.mode_)
      return (mode_ < b.mode_) ? -1 : 1;

    cFlagMapIterator it1 =   flagMap.begin();
    cFlagMapIterator it2 = b.flagMap.begin();

    while(it1 != flagMap.end()){
      const std::string &s1 = it1->second;
      const std::string &s2 = it2->second;

      const int cmp = s1.compare(s2);

      if(cmp)
        return cmp;

      ++it1;
      ++it2;
    }

    return 0;
  }
  //==================================
};
