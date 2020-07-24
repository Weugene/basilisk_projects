
#include "windows.h"
#include "stdio.h"

#include "hf_error.h"

#ifdef WIN32

	void	ReportError(const char* msg, int errorCode, const char* fName, int lineNum) {
		char	buffer[MAX_PATH];
		_snprintf(buffer, MAX_PATH - 1, "%s (0x%x) at %s : %d\n", msg, errorCode, fName, lineNum);
		printf(buffer);
		OutputDebugString(buffer);
	}

#endif