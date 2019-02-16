<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
          "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    $if(refresh)$
    <meta http-equiv="refresh" content="$refresh$">
    $endif$
    $if(feed)$
    <link href="$base$/_feed/" type="application/atom+xml" rel="alternate" title="Sitewide ATOM Feed" />
    <link href="$base$/_feed$pageUrl$" type="application/atom+xml" rel="alternate" title="This page's ATOM Feed" />
    $endif$
    <title>$wikititle$ - $pagetitle$</title>
    $if(printable)$
    <link href="$base$/css/print.css" rel="stylesheet" media="all" type= "text/css" />
    $else$
    <link href="$base$/css/custom.css" rel="stylesheet" media="screen, projection" type="text/css" />
    <link href="$base$/css/print.css" rel="stylesheet" media="print" type= "text/css" />
    $endif$
    <!--[if IE]><link href="$base$/css/ie.css" rel="stylesheet" media="screen, projection" type="text/css" /><![endif]-->
    <link rel="stylesheet" href="$base$/css/codemirror.css"/>
    <link rel="stylesheet" href="$base$/css/basilisk.css"/>
    <link rel="stylesheet" href="$base$/codemirror/addon/display/fullscreen.css"/>
    <link rel="stylesheet" href="$base$/codemirror/addon/fold/foldgutter.css" />
    <script src="$base$/codemirror/lib/codemirror.js"></script>
    <script src="$base$/codemirror/mode/clike/clike.js"></script>
    <script src="$base$/codemirror/mode/markdown/markdown.js"></script>
    <script src="$base$/codemirror/addon/display/fullscreen.js"></script>
    <script src="$base$/codemirror/addon/edit/matchbrackets.js"></script>
    <script src="$base$/codemirror/keymap/emacs.js"></script>
    <script src="$base$/codemirror/addon/comment/comment.js"></script>
    <script src="$base$/codemirror/addon/dialog/dialog.js"></script>
    <script src="$base$/codemirror/addon/search/searchcursor.js"></script>
    <script src="$base$/codemirror/addon/search/search.js"></script>

    <script src="$base$/codemirror/addon/fold/foldcode.js"></script>
    <script src="$base$/codemirror/addon/fold/foldgutter.js"></script>
    <script src="$base$/codemirror/addon/fold/brace-fold.js"></script>
    <script src="$base$/codemirror/addon/fold/xml-fold.js"></script>
    <script src="$base$/codemirror/addon/fold/markdown-fold.js"></script>
    <script src="$base$/codemirror/addon/fold/comment-fold.js"></script>
  </head>
  <body>
    <div id="doc3" class="yui-t1">
        <div id="yui-main">
          <div id="maincol" class="yui-b">
            $userbox()$
            $tabs$ 
            $content()$
            $footer()$
          </div>
        </div>
        <div id="sidebar" class="yui-b first">
          $logo()$
          $if(sitenav)$
            $sitenav()$
          $endif$
          $if(pagetools)$
            $pagetools()$
          $endif$
          $if(markuphelp)$
            $markuphelp()$
          $endif$
        </div>
    </div>
    $javascripts$
    $expire()$
    $getuser()$
  </body>
</html>
