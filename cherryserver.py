"""
"""

import cherrypy

cherrypy.config.update({
    'server.socket_host': '152.14.14.18',
    'server.socket_port': 8080})

class HelloWorld(object):
    def index(self):
        return 'ohai'
    index.exposed = True

if __name__ == '__main__':
    cherrypy.quickstart(HelloWorld())

