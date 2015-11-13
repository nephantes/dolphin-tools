#!/usr/bin/python
import logging
import urllib,urllib2
import json
import os
import re
import cgi
import warnings
import sys
import time

class funcs:
    re_string = re.compile(r'(?P<htmlchars>[<&>])|(?P<space>^[ \t]+)|(?P<lineend>\r\n|\r|\n)|(?P<protocal>(^|\s)((http|ftp)://.*?))(\s|$)', re.S|re.M|re.I)

    def plaintext2html(self, text, tabstop=4):
        def do_sub(m):
            c = m.groupdict()
            if c['htmlchars']:
                return cgi.escape(c['htmlchars'])
            if c['lineend']:
                return '<br>'
            elif c['space']:
                t = m.group().replace('\t', '&nbsp;'*tabstop)
                t = t.replace(' ', '&nbsp;')
                return t
            elif c['space'] == '\t':
                return ' '*tabstop;
            else:
                url = m.group('protocal')
                if url.startswith(' '):
                    prefix = ' '
                    url = url[1:]
                else:
                    prefix = ''
                last = m.groups()[-1]
                if last in ['\n', '\r', '\r\n']:
                    last = '<br>'
                return '%s<a href="%s">%s</a>%s' % (prefix, url, url, last)
        return re.sub(self.re_string, do_sub, text)
    
    def queryAPI(self, url, data, name, logging):
        opener = urllib2.build_opener(urllib2.HTTPHandler())
        trials=0
        while trials<5:
           try:
              mesg = opener.open(url, data=data).read()
              trials=10
           except:
              print "Couldn't connect to dolphin server (%s)"%trials
              logging.info("Couldn't connect to dolphin server (%s)"%trials)
              time.sleep(15)
           trials=trials+1
        ret=str(json.loads(mesg))
        logging.info("%s:%s"%(name,ret))
    
        if (ret.startswith("ERROR")):
              logging.info("%s:%s"%(name,ret))
              print name + ":" + ret + "\n"
              sys.exit(2);
        return ret
        

def main():

   try:
      print "HERE" 
   except:
      sys.exit(2)

if __name__ == "__main__":
    main()
