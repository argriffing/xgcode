"""
Access the CIPRES Science Gateway (CSG).

Try to use its REST API.
This should unsurprisingly act a lot like hpcutil,
because the CSG runs on a supercomputer in San Diego.
Some information is available at
http://www.phylo.org/rest/
Test using
Apparently the CIPRES REST API was never finished,
so neither is this module because there is nothing for it to connect to.
"""


g_test_command_line = """\
curl -F "datafile=@mydna.nex" -F email="user@server.com"
-F "analysis=ML GTR+G" -F tool=RAxML -F datatype=DNA_DATATYPE
-D h.txt http://8ball.sdsc.edu:8888/cipres-web/restapi/job
"""

g_request_body = """\
POST / HTTP/1.1
User-Agent: Jakarta Commons-HttpClient/3.1-rc1
Host:127.0.0.1:3333 Content-Length: 452
Content-Type: multipart/form-data; boundary=frontier

--frontier
Content-Disposition: form-data; name="datafile";
filename="datafile.nex"
Content-Type: application/octent-stream; charset=ISO-8859-1
Content-Transfer-Encoding:binary

[content of datafile here....]
--frontier
Content-Disposition: form-data; name="configfile";
filename="configfile.xml"
Content-Type: application/octent-stream; charset=ISO-8859-1
Content-Transfer-Encoding:binary

[content of configfile here....]
--frontier
Content-Disposition: form-data; name="email"
Content-Type:text/plain; charset=US-ASCII
Content-Transfer-Encoding:8bit

user@server.com
--frontier
Content-Disposition: form-data; name="datatype"
Content-Type:text/plain; charset=US-ASCII
Content-Transfer-Encoding:8bit

DNA_DATATYPE
--frontier
Content-Disposition: form-data; name="analysis"
Content-Type:text/plain; charset=US-ASCII
Content-Transfer-Encoding:8bit

MP
--frontier
Content-Disposition: form-data; name="tool"
Content-Type:text/plain; charset=US-ASCII
Content-Transfer-Encoding:8bit
PAUP
--frontier--
"""

def main():
    pass

if __name__ == '__main__':
    main()
