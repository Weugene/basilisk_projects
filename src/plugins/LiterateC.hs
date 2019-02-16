{-# LANGUAGE ForeignFunctionInterface #-}
module LiterateC (plugin) where

-- This plugin reformats C comments and code blocks prior to parsing by pandoc.

import Foreign
import Foreign.C
import System.IO.Unsafe as Unsafe
import Control.Monad.State (get)
import Network.Gitit.Interface

plugin :: Plugin
plugin = PreParseTransform literate

literate :: String -> PluginM String
literate x = do
  ctx <- get
  mbuser <- askUser
  let username = case mbuser of
                   Nothing  -> "???"
                   Just u   -> uUsername u
  return $ toLiterate x (pgPageName $ ctxLayout ctx) username

-- from http://stackoverflow.com/questions/16888135/how-to-write-a-pure-string-to-string-function-in-haskell-ffi-to-c
-- and also http://book.realworldhaskell.org/read/interfacing-with-c-the-ffi.html

-- the actual transformation is done in C, Haskell is just too hard...
foreign import ccall unsafe "literate"
     literate_c :: CString -> CString -> CString -> IO CString

-- from http://rosettacode.org/wiki/Call_a_foreign-language_function#Haskell
toLiterate :: String -> String -> String -> String
toLiterate s pname username = Unsafe.unsafePerformIO $
      withCString pname $ \cname -> do
      withCString username $ \cuser -> do
      withCString s -- marshall the Haskell string into a C string...
        (\cs -> -- ... and name it cs
         do s2 <- literate_c cs cname cuser
            s2_hs <- peekCString s2 -- marshall the C string called s2 into a Haskell string named s2_hs
	    free s2
	    return s2_hs) -- cs is automatically freed by withCString once done
