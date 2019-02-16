<div id="content">
  $if(revision)$
    <h2 class="revision">Revision $revision$</h2>
  $endif$
  <h1 class="pageTitle"><a href="$pageUrl$">$pagetitle$</a></h1>
  $if(messages)$
    $messages()$
  $endif$
  $content$
  <script type="text/javascript">
    CodeMirror.keyMap.emacs.Enter = "newlineAndIndent";
    CodeMirror.extendMode("clike", {electricChars: "{};"});
    var editor = 
    	CodeMirror.fromTextArea(document.getElementById("editedText"), {
          matchBrackets: true,
          mode: "text/x-basiliskc",
	  theme: "basilisk",
          tabSize: 8,
	  electricChars: true,
	  lineWrapping: true,
          foldGutter: true,
          keyMap: "emacs",
      	  extraKeys: {
            "F11": function(cm) {
              cm.setOption("fullScreen", !cm.getOption("fullScreen"));
            },
	    "Esc": function(cm) {
              if (cm.getOption("fullScreen")) cm.setOption("fullScreen", false);
            },
	    "Ctrl-E": function(cm){ cm.foldCode(cm.getCursor()); }
          },
	  gutters: ["CodeMirror-foldgutter"]
        });
    editor.on("blur", function() {editor.save()});
    update();
    var pending;
    editor.on("change", function() {
	clearTimeout(pending);
	pending = setTimeout(update, 400);
    });
    function looksLikeBasiliskC(code) {
	return /^[/][*][*]$$/.test(code);
    }
    function update() {
	editor.setOption("mode", looksLikeBasiliskC(editor.getValue()) ? 
			 "text/x-basiliskc" : "markdown");
    }
  </script>
</div>
