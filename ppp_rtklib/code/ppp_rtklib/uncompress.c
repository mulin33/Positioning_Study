#include "rtklib.h"

/* execute command -------------------------------------------------------------
* execute command line by operating system shell
* args   : char   *cmd      I   command line
* return : execution status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int execcmd(const char *cmd)
{
#ifdef WIN32
	PROCESS_INFORMATION info;
	STARTUPINFO si = { 0 };
	DWORD stat;
	char cmds[1024];

	trace(3, "execcmd: cmd=%s\n", cmd);

	si.cb = sizeof(si);
	sprintf(cmds, "cmd /c %s", cmd);
	if (!CreateProcess(NULL, (LPTSTR)cmds, NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL,
		NULL, &si, &info)) return -1;
	WaitForSingleObject(info.hProcess, INFINITE);
	if (!GetExitCodeProcess(info.hProcess, &stat)) stat = -1;
	CloseHandle(info.hProcess);
	CloseHandle(info.hThread);
	return (int)stat;
#else
	trace(3, "execcmd: cmd=%s\n", cmd);

	return system(cmd);
#endif
}

/* uncompress file -------------------------------------------------------------
* uncompress (uncompress/unzip/uncompact hatanaka-compression/tar) file
* args   : char   *file     I   input file
*          char   *uncfile  O   uncompressed file
* return : status (-1:error,0:not compressed file,1:uncompress completed)
* note   : creates uncompressed file in tempolary directory
*          gzip and crx2rnx commands have to be installed in commands path
*-----------------------------------------------------------------------------*/
extern int uncompress(const char *file, char *uncfile)
{
	int stat = 0;
	char *p, cmd[2048] = "", tmpfile[1024] = "";
#ifdef WIN32
	char buff[1024], *fname, *dir = "";
#endif

	trace(3, "uncompress: file=%s\n", file);

	strcpy(tmpfile, file);
	if (!(p = strrchr(tmpfile, '.'))) return 0;

	/* uncompress by gzip */
	if (!strcmp(p, ".z") || !strcmp(p, ".Z") ||
		!strcmp(p, ".gz") || !strcmp(p, ".GZ") ||
		!strcmp(p, ".zip") || !strcmp(p, ".ZIP")) {

		strcpy(uncfile, tmpfile); uncfile[p - tmpfile] = '\0';
		sprintf(cmd, "gzip -f -d -c \"%s\" > \"%s\"", tmpfile, uncfile);

		if (execcmd(cmd)) {
			remove(uncfile);
			return -1;
		}
		strcpy(tmpfile, uncfile);
		stat = 1;
	}
	/* extract tar file */
	if ((p = strrchr(tmpfile, '.')) && !strcmp(p, ".tar")) {

		strcpy(uncfile, tmpfile); uncfile[p - tmpfile] = '\0';
#ifdef WIN32
		strcpy(buff, tmpfile);
		fname = buff;
		if ((p = strrchr(buff, '\\'))) {
			*p = '\0'; dir = fname; fname = p + 1;
		}
		sprintf(cmd, "set PATH=%%CD%%;%%PATH%% & cd /D \"%s\" & tar -xf \"%s\"",
			dir, fname);
#else
		sprintf(cmd, "tar -xf \"%s\"", tmpfile);
#endif
		if (execcmd(cmd)) {
			if (stat) remove(tmpfile);
			return -1;
		}
		if (stat) remove(tmpfile);
		stat = 1;
	}
	/* extract hatanaka-compressed file by cnx2rnx */
	else if ((p = strrchr(tmpfile, '.')) && strlen(p)>3 && (*(p + 3) == 'd' || *(p + 3) == 'D')) {

		strcpy(uncfile, tmpfile);
		uncfile[p - tmpfile + 3] = *(p + 3) == 'D' ? 'O' : 'o';
		sprintf(cmd, "crx2rnx < \"%s\" > \"%s\"", tmpfile, uncfile);

		if (execcmd(cmd)) {
			remove(uncfile);
			if (stat) remove(tmpfile);
			return -1;
		}
		if (stat) remove(tmpfile);
		stat = 1;
	}
	trace(3, "uncompress: stat=%d\n", stat);
	return stat;
}