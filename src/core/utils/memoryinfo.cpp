/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2013 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "voreen/core/utils/memoryinfo.h"

#include "voreen/core/utils/stringutils.h"

#include "tgt/logmanager.h"

#ifdef WIN32
#include "windows.h"
#include "psapi.h"
#endif

namespace voreen {

const std::string MemoryInfo::loggerCat_("voreen.MemoryInfo");

uint64_t MemoryInfo::getTotalVirtualMemory() {
#ifdef WIN32
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&memInfo))
        return static_cast<uint64_t>(memInfo.ullTotalPageFile);
    else
        return 0;
#else
    return 0;
#endif
}

uint64_t MemoryInfo::getAvailableVirtualMemory() {
#ifdef WIN32
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&memInfo))
        return static_cast<uint64_t>(memInfo.ullAvailPageFile);
    else
        return 0;
#else
    return 0;
#endif
}

uint64_t MemoryInfo::getVirtualMemoryUsedByCurrentProcess() {
#ifdef WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return static_cast<uint64_t>(pmc.PagefileUsage);
    else
        return 0;
#else
    return 0;
#endif
}

uint64_t MemoryInfo::getTotalPhysicalMemory() {
#ifdef WIN32
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&memInfo))
        return static_cast<uint64_t>(memInfo.ullTotalPhys);
    else
        return 0;
#else
    return 0;
#endif
}

uint64_t MemoryInfo::getAvailablePhysicalMemory() {
#ifdef WIN32
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&memInfo))
        return static_cast<uint64_t>(memInfo.ullAvailPhys);
    else
        return 0;
#else
    return 0;
#endif
}

uint64_t MemoryInfo::getPhysicalMemoryUsedByCurrentProcess() {
#ifdef WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return static_cast<uint64_t>(pmc.WorkingSetSize);
    else
        return 0;
#else
    return 0;
#endif
}

std::string MemoryInfo::getTotalMemoryAsString() {
    return "Total CPU RAM (physical/virtual): " +
            formatMemorySize(getTotalPhysicalMemory()) + " / " +
            formatMemorySize(getTotalVirtualMemory());
}

std::string MemoryInfo::getAvailableMemoryAsString() {
    return "Available CPU RAM (physical/virtual): " +
           formatMemorySize(getAvailablePhysicalMemory()) + " / " +
           formatMemorySize(getAvailableVirtualMemory());
}

std::string MemoryInfo::getProcessMemoryUsageAsString() {
    return "CPU RAM used by process (physical/): " +
           formatMemorySize(getPhysicalMemoryUsedByCurrentProcess()) + " / " +
           formatMemorySize(getVirtualMemoryUsedByCurrentProcess());
}

} // namespace
