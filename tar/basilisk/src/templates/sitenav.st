<div class="sitenav">
  <fieldset>
    <legend>Site</legend>
    <ul>
      <li><a href="$base$/">Front page</a></li>
      <li><a href="$base$/_index">All pages</a></li>
      <li><a href="$base$/_activity">Recent activity</a></li>
      $if(wikiupload)$
        <li><a href="$base$/_upload">Upload a file</a></li>
      $endif$
      $if(feed)$
      <li><a href="$base$/_feed/" type="application/atom+xml" rel="alternate" title="ATOM Feed">Atom feed</a> <img alt="feed icon" src="$base$/img/icons/feed.png"/></li>
      $endif$
      <li><a href="$base$/Help">Help</a></li>
    </ul>
    <form action="/_search" method="post" id="searchform">
     <input type="text" name="patterns" id="patterns"/>
     <input type="submit" name="search" id="search" value="Search"/>
    </form>
  </fieldset>
  <fieldset>
    <legend>Documentation</legend>
    <ul>
      <li><a href="$base$/Tutorial">Tutorial</a></li>
      <li><a href="$base$/src/INSTALL">Installation</a></li>
      <li><a href="$base$/Basilisk C">Basilisk C</a></li>
      <li><a href="$base$/src/README">Solvers and functions</a></li>
      <li><a href="$base$/src/examples/README">Examples</a></li>
      <li><a href="http://groups.google.com/d/forum/basilisk-fr">User forum</a></li>
    </ul>
  </fieldset>
  <fieldset>
    <legend>Development</legend>
    <ul>
      <li><a href="$base$/src/">Source code</a></li>
      <li><a href="https://hub.darcs.net/basilisk/basilisk/browse/src">Darcs Hub</a></li>
      <li><a href="$base$/src/test/README">Test cases</a></li>
      <li><a href="$base$/sandbox/bugs/README">Bug reports</a></li>
      <li><a href="$base$/sandbox/">Play in the sandbox</a></li>
    </ul>
  </fieldset>
</div>
