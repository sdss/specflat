;+
; NAME:
;   specflat_version
; PURPOSE:
;   Return the version name for the specflat product
; CALLING SEQUENCE:
;   vers = specflat_version([/truncate])
; OUTPUTS:
;   vers       - Version name for the product specflat
;   truncate to first word of vers if /truncate keyword set
; COMMENTS:
;   Depends on shell script in $SPECFLAT_DIR/bin
;-
;------------------------------------------------------------------------------
function specflat_version, truncate=truncate
   cmd = "specflat_version"
   if keyword_set(truncate) then cmd = [cmd,'-t']
   spawn, cmd, stdout, /noshell
   specflat_version = stdout[0]
   return, specflat_version
end
;------------------------------------------------------------------------------
