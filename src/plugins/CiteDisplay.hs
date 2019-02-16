module CiteDisplay (plugin) where

{-
* Renders bibliography data inside a code block using citeproc-hs.
  The data can come from two sources:  it can be included inline
  like so

~~~ {.bib}
@article {Wiles95,
    AUTHOR = {Wiles, Andrew},
     TITLE = {Modular elliptic curves and Fermats last theorem},
   JOURNAL = {Ann. of Math. (2)},
  FJOURNAL = {Annals of Mathematics. Second Series},
    VOLUME = {141},
      YEAR = {1995},
    NUMBER = {3},
     PAGES = {443--551},
      ISSN = {0003-486X},
     CODEN = {ANMAAH},
   MRCLASS = {11G05 (11D41 11F11 11F80 11G18)},
  MRNUMBER = {MR1333035 (96d:11071)},
MRREVIEWER = {Karl Rubin},
       DOI = {10.2307/2118559},
       URL = {http://dx.doi.org/10.2307/2118559},
}
~~~

  or the data can be read from a file like so

~~~ {.bib file="filename"}
Wiles95
Erdos90
~~~

  When reading from a file, the block contains the keys to display from the file
  with keys seperated by spaces or newlines.  If no keys are provided the entire
  file is displayed.

* Each bib entry is labeled in HTML by the key.  So you can link to an entry
  from anywhere else in gitit by using [some text](page name#Wiles95)  Currently
  this only works with HTML.

* If the entry contains a url this plugin adds a link.  You can edit the url
  to point to a file inside gitit or to an external site.  Essentially, everything
  in the URL is put into the parens following a link, i.e. the added link is
  [article](text from URL) so any link type supported by pandoc can be stuck
  into the URL of the reference.

* This plugin uses the citeproc-hs package to parse and render the data. Contrary
  to what the citeproc-hs documentation says, you do not need pandoc compiled
  with support for citeproc-hs.

* This plugin supports any imput format supported by citeproc-hs, which is either
  just "mods" or "mods", "bibtex", "biblatex", "ris", "endnote", "endnotexml", "isi",
  "medline", and "copac" depending on how citeproc-hs was compiled.
  By default, the format is "bibtex" but that can be changed like so

  ~~~ { .bib format="mods" }
  some stuff here
  ~~~

* The output format is controlled by a CSL file in the static directory.  See
  http://www.zotero.org/styles for a large number of style files.  By default,
  the file "bibstyle.csl" from the static directory is used, but you can override
  this for each block like so

  ~~~ {.bib style="filename.csl"}
  some stuff here
  ~~~

  where filename is relative to the static directory.


* If you want to add custom processing to the output, edit the transformOutput
  function below.
-}

import Network.Gitit.Interface

import System.Directory (getTemporaryDirectory, removeFile)
import System.IO (openTempFile, hPutStrLn, hClose)
import System.FilePath ((</>))
import Data.Maybe (fromMaybe)
import Control.Monad (liftM)
import Data.Char (toLower)
import Control.Exception
import Data.FileStore (FileStore, retrieve)
import qualified Data.List as L
import qualified Data.Map as M
import qualified Text.CSL as CSL
import qualified Text.ParserCombinators.Parsec as P
import Text.ParserCombinators.Parsec ((<|>))

-------------------------------------------------------------------------------
--- Parsing and rendering bibliography data
-------------------------------------------------------------------------------

defaultFormat = "bibtex"
defaultStyleFile = "bibstyle.csl"

plugin :: Plugin
plugin = mkPageTransformM transformBlock

data RawInput = RawInput String String deriving (Show,Eq) -- first is the format, second is the data.

transformBlock :: Block -> PluginM Block
transformBlock (CodeBlock (_, classes, namevals) contents) | "bib" `elem` classes = do
  cfg <- askConfig
  fstore <- askFileStore
  let format = fromMaybe defaultFormat $ L.lookup "format" namevals
      styleFile = fromMaybe defaultStyleFile $ L.lookup "style" namevals

  refs <- liftIO $ try $ parseBiblio fstore format contents $ L.lookup "file" namevals
  style <- liftIO $ try $ parseCSLFile cfg styleFile

  case refs of
    Left err -> return $ BlockQuote [Para [Str $ "Error parsing bib data: " ++ show err]]
    Right (refs',raw) -> case style of
                                Left err -> return $ BlockQuote [Para [Str $ "Error parsing style file: " ++ show err]]
                                Right style' -> return $ renderBibData raw refs' style'

transformBlock x = return x

--- Given parsed references and a style, render them into a definition list
renderBibData :: RawInput -> [CSL.Reference] -> CSL.Style -> Block
renderBibData raw refs style = DefinitionList $ L.map (\(x,y,_) -> (x,[Plain y])) $ transformOutput raw output
  where formatted = L.concatMap (CSL.processBibliography style . (:[])) refs
        rendered = L.map (read . CSL.renderPandoc style) formatted
        output = zip3 (map ((:[]) . Str . CSL.citeId) refs) rendered refs

--- Parses the bib data.  Returns the list of references and the unparsed data.
parseBiblio :: FileStore -> String -> String -> Maybe String -> IO ([CSL.Reference],RawInput)
parseBiblio fstore format contents Nothing = do
  tdir <- getTemporaryDirectory
  (f,h) <- openTempFile tdir "gitit.biblio"
  let tfile = tdir </> f
  hPutStrLn h contents
  hClose h

  refs <- CSL.readBiblioFile tfile format
  removeFile f
  return (refs, RawInput format contents)

parseBiblio fstore format contents (Just filename) = do
  str <- retrieve fstore filename Nothing
  (refs,raw) <- parseBiblio fstore format str Nothing
  if L.null keys
    then return (refs, raw)
    else return (L.filter (\x -> CSL.citeId x `elem` keys) refs, raw)
 where keys = words contents

parseCSLFile :: Config -> String -> IO CSL.Style
parseCSLFile cfg s = CSL.readCSLFile (staticDir cfg </> s)


-------------------------------------------------------------------------------
--- Custom transforms on the output after citeproc-hs has rendered it
-------------------------------------------------------------------------------

-- The output is a definition list.  The first element of the tuple is the item being defined,
-- by default the citeId.  The second element of the tuple is the display of the reference, and the
-- third element is the parsed Reference from citeproc-hs.
type Output = ([Inline], [Inline], CSL.Reference)

transformOutput :: RawInput -> [Output] -> [Output]
transformOutput rawdata = L.map (addLink . addURL . addMRNumber rawdata)

--- Add a link to the article if the url is not null.
addURL :: Output -> Output
addURL (d,c,r) | (not $ L.null $ CSL.url r) = (d, c ++ [Str " ", Link [Str "article"] (CSL.url r,[])], r)
addURL x = x

--- Add a label to this entry so it can be linked from other pages in gitit.
addLink :: Output -> Output
addLink (d,c,r) = ([RawInline $ "<a name=\"" ++ CSL.citeId r ++ "\"></a>"] ++ d, c, r)

--- Sadly the MR number is not parsed by citeproc-hs so the MR number does not appear in the
--- Reference anywhere.  Rather than hack citeproc-hs I just added a small bibtex parser
--- since the MR number will probably only be included if the bibtex was copied from
--- MathSciNet

addMRNumber :: RawInput -> Output -> Output
addMRNumber (RawInput format rawdata) (d,c,r) | format == "bibtex" =
     case M.lookup (CSL.citeId r, "mrnumber") bibitems of
       Just x -> (d, c ++ mkURL x, r)
       Nothing -> (d, c, r)
  where bibitems = case P.parse bibParser "" rawdata of
                     Left err -> M.empty -- todo: add error reporting?
                     Right x -> x
        mkURL x | length x > 2 = [Str " ", Link [Str "MathSciNet"] (mathSciNet ++ (drop 2 $ head $ words x), [])]
        mkURL _ = []
        mathSciNet = "http://www.ams.org/mathscinet/search/publdoc.html?pg1=MR&s1="
addMRNumber _ x = x

bibParser :: P.Parser (M.Map (String,String) String)
bibParser = do
  x <- P.sepEndBy bibEntry P.spaces
  P.eof
  return $ foldl (\m (k,l) -> M.union m $ M.mapKeys (\k' -> (k,k')) $ M.fromList l) M.empty x

bibEntry :: P.Parser (String, [(String,String)])
bibEntry = do
  P.char '@'
  P.many $ P.noneOf "{"
  P.char '{'
  key <- P.many $ P.noneOf ","
  P.many $ P.oneOf " \t\n,"
  attrs <- P.sepEndBy bibAttr $ P.many1 $ P.oneOf " \t\n,"
  P.char '}'
  return (key, attrs)

bibAttr :: P.Parser (String, String)
bibAttr = do
  key <- P.many (P.letter <|> P.digit)
  P.spaces
  P.char '='
  P.spaces
  P.char '{'
  val <- bibVal
  P.char '}'
  return (map toLower key, val)

bibVal :: P.Parser String
bibVal = liftM concat $ P.many1 (bibValMatched <|> (liftM (:[]) (P.noneOf "{}")))

bibValMatched :: P.Parser String
bibValMatched = P.between (P.char '{') (P.char '}') bibVal
