/******************************************************************************
 * src/wxClickText.h
 *
 * Implements a wxStaticText label which allows a tooltip and click events.
 *
 ******************************************************************************
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef WXCLICKTEXT_H
#define WXCLICKTEXT_H

#include <wx/stattext.h>

class wxClickText : public wxStaticText
{
public:

    wxClickText(wxWindow *parent, wxWindowID id, const wxString &label,
                const wxPoint& pos = wxDefaultPosition,
                const wxSize& size = wxDefaultSize,
                int style = 0, const wxString& name = _T("staticText"));

    virtual ~wxClickText();

    void OnMouseLeftDownEvent(wxMouseEvent& event);

    DECLARE_EVENT_TABLE()
};

#endif // WXCLICKTEXT_H
